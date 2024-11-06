"""
CLI tool to place lipids in membrane domains based on local curvature preferences.
Uses domain_input.txt format for lipid specifications and generates input.str for next steps.

Example domain_input.txt:
; domain lipid percentage c0 density
0 POPC .5 0.179 0.64
2 POPG .5 0.629 0.64
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence

from ..core.point import Point

logger = logging.getLogger(__name__)

@dataclass
class LipidSpec:
    """Specification for a lipid type and its properties"""
    domain_id: int
    name: str
    percentage: float
    curvature: float
    density: float

def parse_lipid_file(file_path: Path) -> List[LipidSpec]:
    """Parse lipid specification file into structured data"""
    lipids = []
    with open(file_path) as f:
        for line in f:
            if line.strip() and not line.startswith(';'):
                domain_id, name, percentage, curvature, density = line.split()
                lipids.append(LipidSpec(
                    domain_id=int(domain_id),
                    name=name,
                    percentage=float(percentage),
                    curvature=float(curvature),
                    density=float(density)
                ))

    total = sum(lipid.percentage for lipid in lipids)
    if not np.isclose(total, 1.0, atol=0.01):
        raise ValueError(f"Lipid percentages must sum to 1.0 (got {total:.2f})")

    return lipids

def write_input_str(lipids: Sequence[LipidSpec], output_file: Path, old_input: Optional[Path] = None) -> None:
    """Write input.str file for TS2CG"""
    # Get existing content from old input if it exists
    existing_content = []
    if old_input and old_input.exists():
        with open(old_input) as f:
            in_lipids = False
            for line in f:
                if "[Lipids List]" in line:
                    in_lipids = True
                elif in_lipids and line.strip().startswith("["):
                    in_lipids = False
                elif not in_lipids and line.strip():
                    existing_content.append(line)

    # Write new file
    with open(output_file, 'w') as f:
        f.write("[Lipids List]\n")
        for lipid in lipids:
            f.write(f"Domain {lipid.domain_id}\n")
            f.write(f"{lipid.name} 1 1 {lipi.density}\n")
            f.write("End\n")
        f.write("\n")
        if existing_content:
            f.writelines(existing_content)

def assign_domains(membrane: Point, lipids: Sequence[LipidSpec], layer: str = "both", k_factor: float = 1.0) -> None:
    """Assign lipids to domains based on curvature preferences"""
    # Determine which layers to process
    layers = [membrane.outer]  # Always process outer layer
    if not membrane.monolayer and layer.lower() in ["both", "inner"]:
        if layer.lower() == "both":
            layers.append(membrane.inner)
        elif layer.lower() == "inner":
            layers = [membrane.inner]

    for membrane_layer in layers:
        layer_name = "outer" if membrane_layer is membrane.outer else "inner"
        logger.info(f"Processing {layer_name} layer")

        n_points = len(membrane_layer.ids)
        curvatures = membrane_layer.mean_curvature

        # Calculate target counts and initialize domains
        target_counts = [int(lipid.percentage * n_points) for lipid in lipids]
        target_counts[-1] += n_points - sum(target_counts)  # Adjust to match total points
        new_domains = np.full(n_points, -1)

        # Mask for tracking unassigned points
        unassigned = np.ones(n_points, dtype=bool)

        # Assign domains
        for i, (lipid, count) in enumerate(zip(lipids, target_counts)):
            if count == 0:
                continue

            # Calculate curvature preference scores for unassigned points
            delta = curvatures[unassigned] - lipid.curvature
            scores = np.exp(-k_factor * delta * delta)

            # Select points with highest scores
            best_indices = np.argpartition(scores, -count)[-count:]
            points_to_assign = np.where(unassigned)[0][best_indices]

            # Assign domain and update mask
            new_domains[points_to_assign] = lipid.domain_id
            unassigned[points_to_assign] = False

        # Update membrane
        membrane_layer.domain_ids = new_domains

        # Log results
        for lipid in lipids:
            actual_count = np.sum(new_domains == lipid.domain_id)
            logger.info(f"{lipid.name}: {actual_count/n_points*100:.1f}% "
                       f"(target: {lipid.percentage*100:.1f}%)")

def DOP(args: List[str]) -> None:
    """Main entry point for Domain Placer tool"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--point-dir', type=Path, default="point",
                       help='Path to membrane point directory (default: point/)')
    parser.add_argument('-i', '--input', type=Path, default="domain_input.txt",
                       help='Path to lipid specification file (default: domain_input.txt)')
    parser.add_argument('-l', '--layer', choices=['both', 'inner', 'outer'],
                       default='both', help='Which membrane layer(s) to modify (default: both)')
    parser.add_argument('-k', '--k-factor', type=float, default=1.0,
                       help='Scaling factor for curvature preference strength (default: 1.0)')
    parser.add_argument('-o', '--output', type=Path,
                       help='Output directory (defaults to input directory)')
    parser.add_argument('-ni', '--new-input', type=Path, default="input.str",
                       help='Path for output input.str file (default: input.str)')
    parser.add_argument('-oi', '--old-input', type=Path,
                       help='Path to existing input.str to preserve additional sections')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')

    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)

    if args.seed is not None:
        np.random.seed(args.seed)

    try:
        membrane = Point(args.point_dir)
        lipids = parse_lipid_file(args.input)

        assign_domains(membrane, lipids, args.layer, args.k_factor)

        write_input_str(lipids, args.new_input, args.old_input)
        logger.info(f"Created input file: {args.new_input}")

        output_dir = args.output or args.point_dir
        membrane.save(output_dir)
        logger.info(f"Updated membrane domains in {output_dir}")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
