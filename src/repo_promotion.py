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
    return json.loads(path.read_text(encoding="utf-8"))


def _iter_manifest_entries(manifest: Dict[str, Any], *, entry_key: str, category: str) -> Iterable[tuple[str, str]]:
    for rel_path in manifest.get(entry_key, []):
        yield category, str(rel_path)


def _resolve_repo_file(root: Path, rel_path: str) -> Path:
    relative = Path(rel_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Manifest path must be relative and confined to the repository: {rel_path}")

    resolved_root = root.resolve()
    candidate = resolved_root / relative
    try:
        candidate.resolve(strict=False).relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"Manifest path resolves outside the repository: {rel_path}") from exc
    return candidate


def build_promotion_plan(
    *,
    source_root: Path,
    target_root: Path,
    manifest: Dict[str, Any],
    assets_manifest: Dict[str, Any] | None = None,
    include_assets: bool = False,
) -> List[PromotionAction]:
    plan: List[PromotionAction] = []
    entries = list(_iter_manifest_entries(manifest, entry_key="sync_files", category="shared_code"))
    if include_assets and assets_manifest is not None:
        entries.extend(_iter_manifest_entries(assets_manifest, entry_key="asset_files", category="demo_asset"))

    removal_paths = [str(path) for path in manifest.get("remove_files", [])]
    overlap = {rel_path for _, rel_path in entries}.intersection(removal_paths)
    if overlap:
        raise ValueError(f"Manifest paths cannot be synchronized and removed together: {sorted(overlap)}")

    for category, rel_path in entries:
        source_path = _resolve_repo_file(source_root, rel_path)
        target_path = _resolve_repo_file(target_root, rel_path)
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

    for rel_path in removal_paths:
        source_path = _resolve_repo_file(source_root, rel_path)
        target_path = _resolve_repo_file(target_root, rel_path)
        if not target_path.exists():
            action = "skip"
            reason = "Obsolete target file is already absent."
        elif not target_path.is_file():
            action = "blocked_remove"
            reason = "Removal entries may target files only; directories are never removed."
        else:
            action = "remove"
            reason = "Target file is explicitly listed as obsolete in the shared manifest."

        plan.append(
            PromotionAction(
                rel_path=rel_path,
                source_path=source_path,
                target_path=target_path,
                category="obsolete_file",
                action=action,
                reason=reason,
            )
        )
    return plan


def apply_promotion_plan(plan: Iterable[PromotionAction]) -> Dict[str, int]:
    counts = {"create": 0, "update": 0, "remove": 0, "skip": 0, "missing_source": 0, "blocked_remove": 0}
    for item in plan:
        counts[item.action] = counts.get(item.action, 0) + 1
        if item.action == "remove":
            if not item.target_path.is_file():
                raise ValueError(f"Refusing to remove a non-file target: {item.target_path}")
            item.target_path.unlink()
            continue
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
    parser = argparse.ArgumentParser(
        description="Promote approved OMIX Exosome Rejuvenation files into the clean public repository."
    )
    parser.add_argument(
        "--target-root",
        type=Path,
        default=None,
        help="Path to the clean public repository. Defaults to the sibling SRSC directory.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Shared code/docs manifest JSON. Defaults to promotion_manifest.json in the source repo.",
    )
    parser.add_argument(
        "--assets-manifest",
        type=Path,
        default=None,
        help="Optional release/demo assets manifest JSON. Defaults to promotion_assets_manifest.json in the source repo.",
    )
    parser.add_argument(
        "--include-assets",
        action="store_true",
        help="Also promote optional release/demo assets listed in the assets manifest.",
    )
    parser.add_argument(
        "--include-demo-data",
        action="store_true",
        help=argparse.SUPPRESS,
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
    assets_manifest_path = (
        args.assets_manifest.resolve()
        if args.assets_manifest is not None
        else source_root / "promotion_assets_manifest.json"
    )
    include_assets = args.include_assets or args.include_demo_data

    manifest = load_promotion_manifest(manifest_path)
    assets_manifest = load_promotion_manifest(assets_manifest_path) if include_assets else None
    plan = build_promotion_plan(
        source_root=source_root,
        target_root=target_root,
        manifest=manifest,
        assets_manifest=assets_manifest,
        include_assets=include_assets,
    )

    logger.info("Source root: %s", source_root)
    logger.info("Target root: %s", target_root)
    logger.info("Shared manifest: %s", manifest_path)
    if include_assets:
        logger.info("Assets manifest: %s", assets_manifest_path)
    for row in _plan_rows(plan):
        logger.info("%s | %s | %s | %s", row["action"].upper(), row["category"], row["rel_path"], row["reason"])

    if args.apply:
        counts = apply_promotion_plan(plan)
        logger.info("Promotion applied with counts: %s", counts)
    else:
        logger.info("Dry run only. Re-run with --apply to synchronize approved files with the clean repo.")


if __name__ == "__main__":
    main()
