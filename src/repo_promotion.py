from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

from .logging_utils import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class PromotionAction:
    rel_path: str
    source_path: Path
    target_path: Path
    category: str
    action: str
    reason: str


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_promotion_manifest(path: Path) -> Dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest.setdefault("sync_files", [])
    manifest.setdefault("optional_demo_files", [])
    return manifest


def _iter_manifest_entries(manifest: Dict[str, Any], include_demo_data: bool) -> Iterable[tuple[str, str]]:
    for rel_path in manifest.get("sync_files", []):
        yield "sync", str(rel_path)
    if include_demo_data:
        for rel_path in manifest.get("optional_demo_files", []):
            yield "demo_data", str(rel_path)


def build_promotion_plan(
    *,
    source_root: Path,
    target_root: Path,
    manifest: Dict[str, Any],
    include_demo_data: bool = False,
) -> List[PromotionAction]:
    plan: List[PromotionAction] = []
    for category, rel_path in _iter_manifest_entries(manifest, include_demo_data):
        source_path = source_root / rel_path
        target_path = target_root / rel_path
        if not source_path.exists():
            plan.append(
                PromotionAction(
                    rel_path=rel_path,
                    source_path=source_path,
                    target_path=target_path,
                    category=category,
                    action="missing_source",
                    reason="Source file is not present in the work repo.",
                )
            )
            continue

        if not target_path.exists():
            plan.append(
                PromotionAction(
                    rel_path=rel_path,
                    source_path=source_path,
                    target_path=target_path,
                    category=category,
                    action="create",
                    reason="Target file does not exist in the clean repo.",
                )
            )
            continue

        if _file_digest(source_path) == _file_digest(target_path):
            action = "skip"
            reason = "Source and target files are identical."
        else:
            action = "update"
            reason = "Source and target files differ."

        plan.append(
            PromotionAction(
                rel_path=rel_path,
                source_path=source_path,
                target_path=target_path,
                category=category,
                action=action,
                reason=reason,
            )
        )
    return plan


def apply_promotion_plan(plan: Iterable[PromotionAction]) -> Dict[str, int]:
    counts = {"create": 0, "update": 0, "skip": 0, "missing_source": 0}
    for item in plan:
        counts[item.action] = counts.get(item.action, 0) + 1
        if item.action not in {"create", "update"}:
            continue
        item.target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item.source_path, item.target_path)
    return counts


def _plan_rows(plan: Iterable[PromotionAction]) -> List[Dict[str, Any]]:
    return [
        {
            **asdict(item),
            "source_path": str(item.source_path),
            "target_path": str(item.target_path),
        }
        for item in plan
    ]


def _default_target_root(source_root: Path) -> Path:
    return source_root.parent / "SRSC"


def main() -> None:
    parser = argparse.ArgumentParser(description="Promote approved SRSC-work files into the clean SRSC repo.")
    parser.add_argument(
        "--target-root",
        type=Path,
        default=None,
        help="Path to the clean SRSC repository. Defaults to the sibling SRSC directory.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Promotion manifest JSON. Defaults to promotion_manifest.json in the source repo.",
    )
    parser.add_argument(
        "--include-demo-data",
        action="store_true",
        help="Also promote demo/sample files listed in optional_demo_files.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply the promotion plan. Without this flag, the script only prints a dry-run plan.",
    )
    args = parser.parse_args()

    source_root = Path(__file__).resolve().parents[1]
    target_root = args.target_root.resolve() if args.target_root is not None else _default_target_root(source_root)
    manifest_path = args.manifest.resolve() if args.manifest is not None else source_root / "promotion_manifest.json"

    manifest = load_promotion_manifest(manifest_path)
    plan = build_promotion_plan(
        source_root=source_root,
        target_root=target_root,
        manifest=manifest,
        include_demo_data=args.include_demo_data,
    )

    logger.info("Source root: %s", source_root)
    logger.info("Target root: %s", target_root)
    logger.info("Manifest: %s", manifest_path)
    for row in _plan_rows(plan):
        logger.info("%s | %s | %s", row["action"].upper(), row["rel_path"], row["reason"])

    if args.apply:
        counts = apply_promotion_plan(plan)
        logger.info("Promotion applied with counts: %s", counts)
    else:
        logger.info("Dry run only. Re-run with --apply to copy approved files into the clean repo.")


if __name__ == "__main__":
    main()
