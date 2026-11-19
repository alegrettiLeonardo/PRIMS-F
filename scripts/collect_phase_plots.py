"""Generate phase comparison figures once immutable terminal analyses exist."""
import argparse
import hashlib
import json
from pathlib import Path

from plot_phase_comparison import plot


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True)
    evidence=[]
    for path in sorted(args.analysis_root.rglob("*.json")):
        try:
            report=json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        if not isinstance(report,dict) or report.get("schema")!="PRIMS-F.phase_comparison.v1":
            continue
        relative=path.relative_to(args.analysis_root)
        target=args.output/("_".join(relative.with_suffix("").parts)+".png")
        if target.exists():
            raise FileExistsError(target)
        plot(report,target)
        evidence.append({"source":str(path),"source_sha256":sha(path),
                         "figure":str(target),"figure_sha256":sha(target)})
    with (args.output/"index.json").open("x") as stream:
        json.dump(evidence,stream,indent=2)
        stream.write("\n")


if __name__=="__main__":
    main()
