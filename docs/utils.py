import os
import re
from logging import info, warning
from pathlib import Path
from shutil import copytree, rmtree
from tempfile import TemporaryDirectory
from typing import Any, Dict, ForwardRef, List, Union

from git import Repo

HERE = Path(__file__).parent


def _fetch_notebooks(repo_url: str) -> None:
    def copy_files(repo_path: Union[str, Path]) -> None:
        repo_path = Path(repo_path)

        for dirname in ["tutorials", "auto_examples", "gen_modules"]:
            rmtree(dirname, ignore_errors=True)  # locally re-cloning
            copytree(repo_path / "docs" / "source" / dirname, dirname)

    def fetch_remote(repo_url: str) -> None:
        info(f"Fetching notebooks from repo `{repo_url}`")
        with TemporaryDirectory() as repo_dir:
            ref = "main"
            repo = Repo.clone_from(repo_url, repo_dir, depth=1, branch=ref)
            repo.git.checkout(ref, force=True)

            copy_files(repo_dir)

    def fetch_local(repo_path: Union[str, Path]) -> None:
        info(f"Fetching notebooks from local path `{repo_path}`")
        repo_path = Path(repo_path)
        if not repo_path.is_dir():
            raise OSError(f"Path `{repo_path}` is not a directory.")

        copy_files(repo_path)

    notebooks_local_path = Path(
        os.environ.get(
            "SQUIDPY_NOTEBOOKS_PATH",
            HERE.absolute().parent.parent.parent / "trackc_notebooks",
        )
    )
    try:
        fetch_local(notebooks_local_path)
    except Exception as e:  # noqa: BLE001
        warning(
            f"Unable to fetch notebooks locally from `{notebooks_local_path}`, reason: `{e}`. Trying remote"
        )
        download = int(os.environ.get("TRACKC_DOWNLOAD_NOTEBOOKS", 1))
        if not download:
            # use possibly old files, otherwise, bunch of warnings will be shown
            info(
                f"Not fetching notebooks from remove because `TRACKC_DOWNLOAD_NOTEBOOKS={download}`"
            )
            return

        fetch_remote(repo_url)
