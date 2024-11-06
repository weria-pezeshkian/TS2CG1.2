
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
from typing import Dict, List, Optional, Sequence
import logging
from dataclasses import dataclass

from .point import Point

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
    specs = []

    try:
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(';'):
                    continue

                parts = line.split()
                spec = LipidSpec(
                    domain_id=int(parts[0]),
                    name=parts[1],
                    percentage=float(parts[2]),
                    curvature=float(parts[3]),
                    density=float(parts[4])
                )
                specs.append(spec)

        total_percentage = sum(spec.percentage for spec in specs)
        if not np.isclose(total_percentage, 1.0, atol=0.01):
            raise ValueError(f"Lipid percentages must sum to 1.0 (got {total_percentage:.2f})")

        return specs

    except (IndexError, ValueError) as e:
        logger.error(f"Error parsing lipid file: {e}")
        raise

def write_input_str(specs: Sequence[LipidSpec], output_file: Path, old_input: Optional[Path] = None) -> None:
    """
    Write input.str file for TS2CG, preserving any non-lipid sections from template file.
    """
    existing_content = []
    if old_input and old_input.exists():
        with open(old_input) as f:
            in_lipids_block = False
            for line in f:
                if "[Lipids List]" in line:
                    in_lipids_block = True
                elif in_lipids_block and line.strip().startswith("["):
                    in_lipids_block = False
                elif not in_lipids_block and line.strip():
                    existing_content.append(line)

    with open(output_file, 'w') as f:
        f.write("[Lipids List]\n")
        for spec in specs:
            f.write(f"Domain {spec.domain_id}\n")
            f.write(f"{spec.name} 1 1 {spec.density}\n")
            f.write("End\n")
        f.write("\n")

        if existing_content:
            f.writelines(existing_content)

def calculate_curvature_scores(curvatures: np.ndarray, specs: Sequence[LipidSpec], k_factor: float) -> np.ndarray:
    """Calculate Boltzmann-weighted scores for each point-lipid combination."""
    n_points = len(curvatures)
    n_specs = len(specs)
    scores = np.zeros((n_points, n_specs))

    for i, spec in enumerate(specs):
        delta = curvatures - spec.curvature  # Remove unnecessary dimension
        scores[:, i] = np.exp(-k_factor * delta * delta)

    return scores

def assign_domains_by_curvature(membrane: Point, specs: Sequence[LipidSpec],
                              layer: str = "both", k_factor: float = 1.0) -> None:
    """
    Assign lipids to domains based on local curvature preferences while strictly
    maintaining specified lipid ratios.
    """
    if membrane.monolayer and layer.lower() == "inner":
        raise ValueError("Cannot modify inner layer - this is a monolayer system")

    # Always process outer membrane
    layers = [membrane.outer]

    # Add inner membrane only if it exists and requested
    if not membrane.monolayer and layer.lower() in ["both", "inner"]:
        if layer.lower() == "both":
            layers.append(membrane.inner)
        else:  # inner only
            layers = [membrane.inner]

    logger.info(f"Processing {'monolayer' if membrane.monolayer else 'bilayer'} system")

    for membrane_layer in layers:
        layer_name = "outer" if membrane_layer is membrane.outer else "inner"
        logger.info(f"Assigning domains for {layer_name} layer")

        n_points = len(membrane_layer.ids)
        curvatures = membrane_layer.mean_curvature

        # Calculate target counts for each lipid type
        target_counts = np.array([int(spec.percentage * n_points) for spec in specs])

        # Adjust last count to ensure total matches n_points
        target_counts[-1] += n_points - target_counts.sum()

        # Calculate initial curvature preference scores for all points
        scores = calculate_curvature_scores(curvatures, specs, k_factor)

        # Initialize domain assignments
        new_domains = np.full(n_points, -1)
        unassigned_mask = np.ones(n_points, dtype=bool)  # Track unassigned points

        # Assign domains while maintaining exact ratios
        for domain_idx, count in enumerate(target_counts):
            if count == 0:
                continue

            # Get scores for unassigned points for this lipid type
            domain_scores = scores[unassigned_mask, domain_idx]

            # Select the best-matching points for this domain
            best_indices = np.argpartition(domain_scores, -count)[-count:]

            # Get the actual point indices from the unassigned points
            unassigned_indices = np.where(unassigned_mask)[0]
            points_to_assign = unassigned_indices[best_indices]

            # Assign domain
            new_domains[points_to_assign] = specs[domain_idx].domain_id

            # Update unassigned mask
            unassigned_mask[points_to_assign] = False

        # Verify all points were assigned
        assert not np.any(new_domains == -1), "Some points were left unassigned"

        # Update membrane layer
        membrane_layer.domain_ids = new_domains

        # Log the actual percentages achieved
        logger.info(f"Achieved lipid percentages in {layer_name} layer:")
        unique, counts = np.unique(new_domains, return_counts=True)
        for spec in specs:
            actual_percent = counts[unique == spec.domain_id][0] / n_points * 100
            logger.info(f"  {spec.name}: {actual_percent:.1f}% (target: {spec.percentage*100:.1f}%)")

def DOP(args: List[str]) -> None:
    """Main entry point for the Domain Placer tool."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-p', '--point-dir', type=Path, default=Path("point"),
                       help='Path to the membrane point directory (default: point/)')
    parser.add_argument('-i', '--input', type=Path, default=Path("domain_input.txt"),
                       help='Path to lipid specification file (default: domain_input.txt)')
    parser.add_argument('-l', '--layer', choices=['both', 'inner', 'outer'],
                       default='both', help='Which membrane layer(s) to modify (default: both)')
    parser.add_argument('-k', '--k-factor', type=float, default=1.0,
                       help='Scaling factor for curvature preference strength (default: 1.0)')
    parser.add_argument('-o', '--output', type=Path,
                       help='Output directory (defaults to modifying input)')
    parser.add_argument('-ni', '--new-input', type=Path, default=Path("input.str"),
                       help='Path for output input.str file (default: input.str)')
    parser.add_argument('-oi', '--old-input', type=Path,
                       help='Path to existing input.str to preserve additional sections')
    parser.add_argument('--seed', type=int,
                       help='Random seed for reproducibility')

    args = parser.parse_args(args)

    logging.basicConfig(level=logging.INFO)

    if args.seed is not None:
        np.random.seed(args.seed)

    try:
        membrane = Point(args.point_dir)
        specs = parse_lipid_file(args.input)

        logger.info("Processing lipids:")
        for spec in specs:
            logger.info(f"  {spec.name}: {spec.percentage*100:.1f}%, "
                       f"curvature={spec.curvature:.3f}, "
                       f"density={spec.density:.3f}")

        assign_domains_by_curvature(
            membrane, specs,
            layer=args.layer,
            k_factor=args.k_factor
        )

        write_input_str(specs, args.new_input, args.old_input)
        logger.info(f"Created input file at {args.new_input}")

        output_dir = args.output or args.point_dir
        membrane.save(output_dir)
        logger.info(f"Successfully updated membrane domains in {output_dir}")

    except Exception as e:
        logger.error(f"Failed to process membrane: {e}", exc_info=True)
        raise
