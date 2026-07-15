from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, List

import yaml
from importlib.resources import as_file, files

try:
    from platformdirs import user_cache_dir, user_data_dir
except ImportError:
    def user_cache_dir(appname: str) -> str:
        return str(Path.home() / ".cache" / appname)

    def user_data_dir(appname: str) -> str:
        return str(Path.home() / ".local" / "share" / appname)


@dataclass(frozen=True)
class ResourceBundle:
    root_dir: Path
    manifest: dict
    resource_files: Dict[str, Path]  # semantic_key -> file path


class ResourceManager:
    """
    - Resources: stored under <resource_base>/resources/<bundle>/v<bundle_version>/
      OR use user_dir directly if provided.
    - Config assets (jquery etc): cached under <cache_base>/config/
      (no versioning, as requested).
    """

    def __init__(
        self,
        *,
        cache_dir: str | None = None,
        resource_root: str | None = None,
        appname: str = "pyfuse",
    ):
        if resource_root:
            self.resource_base = Path(resource_root).expanduser().resolve()
        else:
            self.resource_base = Path(user_data_dir(appname))

        if cache_dir:
            self.cache_base = Path(cache_dir).expanduser().resolve()
        else:
            self.cache_base = Path(user_cache_dir(appname))

        self.resource_base.mkdir(parents=True, exist_ok=True)
        self.cache_base.mkdir(parents=True, exist_ok=True)

        self.resources_cache_base = self.resource_base / "resources"
        self.resources_cache_base.mkdir(parents=True, exist_ok=True)

        self.config_cache_dir = self.cache_base / "config"
        self.config_cache_dir.mkdir(parents=True, exist_ok=True)

    # ---------------- Resources (hg37/hg38/user) ----------------

    def list_cached_versions(self, bundle_name: str) -> List[Path]:
        """Return sorted cached version directories for a bundle."""
        bundle_root = self.resources_cache_base / bundle_name
        if not bundle_root.is_dir():
            return []
        return sorted([p for p in bundle_root.iterdir() if p.is_dir() and p.name.startswith("v")])

    def resolve_from_cache(self, *, default_bundle: str = "default_grch37") -> ResourceBundle:
        """
        Resolve the latest cached bundle version without using packaged resources.
        Raises FileNotFoundError if no cached bundle exists.
        """
        versions = self.list_cached_versions(default_bundle)
        if not versions:
            raise FileNotFoundError(
                f"No cached versions found for bundle '{default_bundle}' under {self.resources_cache_base}"
            )

        root = versions[-1]
        manifest = self._load_manifest_from_dir(root)
        resource_files = self._build_resource_files(root, manifest)
        self._validate(resource_files, manifest)
        return ResourceBundle(root_dir=root, manifest=manifest, resource_files=resource_files)

    def resolve(
        self,
        *,
        user_dir: Optional[str] = None,
                default_bundle: str = "default_grch37",
    ) -> ResourceBundle:
        """
        Returns a ResourceBundle with real filesystem paths.

        Precedence:
          1) user_dir provided -> use that directory (no caching).
                    2) else -> use latest cached/installed default bundle from resource_root.

        Assumes resources/<bundle>/manifest.yaml contains:
          - bundle_version: <int/str>
          - resource_files: {semantic_key: filename, ...}   (required + optional)
          - required_keys: [semantic_key, ...]              (subset of resource_files keys)
        """
        # 1) user-provided directory
        if user_dir:
            root = Path(user_dir).expanduser().resolve()
            manifest = self._load_manifest_from_dir(root)
            resource_files = self._build_resource_files(root, manifest)
            self._validate(resource_files, manifest)
            return ResourceBundle(root_dir=root, manifest=manifest, resource_files=resource_files)

        # 2) managed default bundles from resource_root only (no packaged fallback)
        return self.resolve_from_cache(default_bundle=default_bundle)

    # ---------------- Config assets (packaged) ----------------

    def resolve_config_assets(self, *, force_refresh: bool = False) -> Dict[str, Path]:
        """
        Cache the packaged config directory into <cache_base>/config/
        and return {semantic_key: absolute_path} for config assets defined in
        package config/settings.yaml:

          config_assets: {jquery: jquery-3.6.0.min.js, ...}
          required_config_assets: [jquery, ...]

        No versioning is used for config cache (as requested).
        """
        settings = self._load_packaged_settings()

        mapping = settings.get("config_assets", {})
        required = settings.get("required_config_assets", [])

        if not isinstance(mapping, dict) or not mapping:
            return {}

        if not isinstance(required, list):
            raise ValueError("'required_config_assets' must be a list in settings.yaml")

        if force_refresh and self.config_cache_dir.exists():
            shutil.rmtree(self.config_cache_dir, ignore_errors=True)
            self.config_cache_dir.mkdir(parents=True, exist_ok=True)

        # Materialize config folder once
        if not any(self.config_cache_dir.iterdir()):
            self._materialize_packaged_config_dir(self.config_cache_dir)

        out = {str(k): (self.config_cache_dir / str(fname)) for k, fname in mapping.items()}

        missing = {k: str(out.get(k)) for k in required if k not in out or not out[k].is_file()}
        if missing:
            raise FileNotFoundError(f"Missing required config asset files: {missing}")

        # drop missing optional assets
        out = {k: p for k, p in out.items() if p.is_file()}

        return out

    # ---------------- Internals ----------------

    def _load_packaged_settings(self) -> dict:
        settings_file = files("pyfuse.config") / "settings.yaml"
        try:
            txt = settings_file.read_text(encoding="utf-8")
        except FileNotFoundError as e:
            raise RuntimeError("Packaged settings.yaml not found in package config directory") from e
        return yaml.safe_load(txt) or {}

    def _load_manifest_from_dir(self, root: Path) -> dict:
        mf = root / "manifest.yaml"
        if not mf.is_file():
            raise FileNotFoundError(f"manifest.yaml not found in resource directory: {root}")
        return yaml.safe_load(mf.read_text(encoding="utf-8")) or {}

    def _materialize_packaged_config_dir(self, dest: Path) -> None:
        """
        Copies packaged config assets into dest.
        """
        src = files("pyfuse.config")
        self._copy_traversable_dir_to_dest(src, dest)

    def _copy_traversable_dir_to_dest(self, src_traversable, dest: Path) -> None:
        """
        Copies a Traversable directory (from importlib.resources.files) to dest using
        a staging dir + rename.
        """
        staging = dest.with_name(dest.name + ".staging")

        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True, exist_ok=True)

        with as_file(src_traversable) as src_path:
            shutil.copytree(src_path, staging, dirs_exist_ok=True)

        if dest.exists():
            shutil.rmtree(dest, ignore_errors=True)
        staging.rename(dest)

    def _build_resource_files(self, root: Path, manifest: dict) -> Dict[str, Path]:
        mapping = manifest.get("resource_files")
        if not isinstance(mapping, dict) or not mapping:
            raise ValueError("manifest.yaml must contain a non-empty mapping 'resource_files'")
        return {str(k): (root / str(v)) for k, v in mapping.items()}

    def _validate(self, resource_files: Dict[str, Path], manifest: dict) -> None:
        required_keys = manifest.get("required_keys", [])
        if not isinstance(required_keys, list):
            raise ValueError("manifest.yaml 'required_keys' must be a list")

        missing_keys = [k for k in required_keys if k not in resource_files]
        if missing_keys:
            raise ValueError(f"Required keys missing from manifest 'resource_files': {missing_keys}")

        missing_files = {k: str(resource_files[k]) for k in required_keys if not resource_files[k].is_file()}
        if missing_files:
            raise FileNotFoundError(f"Missing required resource files: {missing_files}")
