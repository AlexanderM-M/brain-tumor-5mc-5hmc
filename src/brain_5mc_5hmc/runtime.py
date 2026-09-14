"""Explicit workspace configuration, stage loading and native build support."""

from pathlib import Path
from types import SimpleNamespace
import importlib
import json
import os
import subprocess

config = SimpleNamespace(
    workspace=None, reference=None, samtools="samtools", workers=4, allow_annotation_download=False
)
_states = {}
STAGES = json.loads(Path(__file__).with_name("stages.json").read_text())


def configure(path=None, *, workspace=None, reference=None, workers=None):
    """Load JSON configuration. Relative paths are relative to the config file."""
    if _states:
        raise RuntimeError("Start a fresh process before changing an initialized workspace.")
    values = {}
    base = Path.cwd()
    if path:
        path = Path(path).resolve()
        base = path.parent
        values = json.loads(path.read_text())
    if workspace is not None:
        values["workspace"] = str(workspace)
    if reference is not None:
        values["reference"] = str(reference)
    if workers is not None:
        values["workers"] = workers
    for key in ("workspace", "reference"):
        value = values.get(key)
        if value:
            value = Path(value).expanduser()
            values[key] = str((base / value).resolve())
    for key in ("workspace", "reference", "samtools", "workers", "allow_annotation_download"):
        if key in values:
            setattr(config, key, values[key])
    if config.workers < 1:
        raise ValueError("workers must be positive")
    return config


def initialized(stage):
    return _states.get(stage) == "ready"


def begin(stage):
    if not config.workspace:
        raise ValueError("Configure a private workspace before initializing analysis stages.")
    if _states.get(stage) == "loading":
        raise RuntimeError("Circular initialization or previous failure: " + stage)
    _states[stage] = "loading"


def finish(stage):
    _states[stage] = "ready"


def initialize(stage):
    spec = STAGES[stage]
    module = importlib.import_module("." + spec["module"], __package__)
    getattr(module, "initialize_" + spec["prefix"])()
    return module


def run(stage):
    spec = STAGES[stage]
    module = importlib.import_module("." + spec["module"], __package__)
    getattr(module, "run_" + spec["prefix"])()


def native_path(name):
    if not config.workspace:
        raise ValueError("Configure a workspace before building or loading native libraries.")
    path = Path(config.workspace) / ".build" / name
    if not path.is_file():
        raise FileNotFoundError(str(path) + "; run brain-5mc-5hmc build-native first")
    return path


def native_source(name):
    return Path(__file__).with_name("native") / name


def build_native():
    """Compile the two C kernels using pysam's bundled htslib headers (Linux)."""
    import pysam

    if not config.workspace:
        raise ValueError("Configure a workspace first")
    out = Path(config.workspace) / ".build"
    out.mkdir(parents=True, exist_ok=True)
    compiler = os.environ.get("CC", "cc")
    for name in ("strengthen_null", "regional_counts"):
        command = [compiler, "-O3", "-std=c99", "-D_GNU_SOURCE", "-fPIC", "-shared"]
        if name == "regional_counts":
            command += ["-I" + directory for directory in pysam.get_include()]
        command += [str(native_source(name + ".c")), "-o", str(out / (name + ".so")), "-lm"]
        subprocess.run(command, check=True)
    return out


def prepare_workspace():
    """Create empty directories and relative aliases; never replace existing files."""
    if not config.workspace:
        raise ValueError("Configure a workspace first")
    root = Path(config.workspace)
    root.mkdir(parents=True, exist_ok=True)
    layout = json.loads(Path(__file__).with_name("layout.json").read_text())
    for alias, target in layout.items():
        src, dst = root / alias, root / target
        src.parent.mkdir(parents=True, exist_ok=True)
        if Path(target).suffix:
            dst.parent.mkdir(parents=True, exist_ok=True)
        else:
            dst.mkdir(parents=True, exist_ok=True)
        if not src.exists() and not src.is_symlink():
            src.symlink_to(os.path.relpath(dst, src.parent))
    for directory in [
        "figures/main/manuscript",
        "supplementary/figures/manuscript",
        "paper/analysis_reports",
    ]:
        (root / directory).mkdir(parents=True, exist_ok=True)
    for stage in ["pooled_analysis", "story_analysis", "focused_followup", "strengthening"]:
        for directory in ["figures", "plans", "results", "validation", "logs", "targeted_bams"]:
            (root / ".analysis" / stage / directory).mkdir(parents=True, exist_ok=True)
    return root


def process_pool(max_workers):
    """Use Linux fork explicitly so initialized stage arrays are inherited."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context

    return ProcessPoolExecutor(max_workers=max_workers, mp_context=get_context("fork"))
