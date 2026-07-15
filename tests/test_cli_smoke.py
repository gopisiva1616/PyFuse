import sys

import pytest

from pyfuse.cli import main


def test_pyfuse_top_level_help(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["pyfuse", "-h"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 0
