"""Command-line entry point; each stage runs in a fresh Python process."""

import argparse
import sys
from . import runtime


def main():
    parser = argparse.ArgumentParser(description="Brain tumour 5mC/5hmC analysis")
    parser.add_argument("--config", help="JSON workspace configuration")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="List the analysis stages")
    sub.add_parser("init-workspace", help="Create an empty private workspace")
    sub.add_parser("build-native", help="Compile the C kernels")
    stage = sub.add_parser("run", help="Run one analysis stage")
    stage.add_argument("stage", choices=list(runtime.STAGES))
    stage.add_argument("stage_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command == "list":
        for name, spec in runtime.STAGES.items():
            print(f"{name:32s} {spec['module']}.{spec['prefix']}")
        return
    if not args.config:
        parser.error("--config is required for this command")
    runtime.configure(args.config)
    if args.command == "init-workspace":
        print(runtime.prepare_workspace())
    elif args.command == "build-native":
        print(runtime.build_native())
    else:
        sys.argv = [args.stage] + args.stage_args
        runtime.run(args.stage)


if __name__ == "__main__":
    main()
