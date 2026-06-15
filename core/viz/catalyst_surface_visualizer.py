#!/usr/bin/env python
"""Catalyst Surface Visualizer with OVITO high-quality rendering only."""

import numpy as np
import argparse
from pathlib import Path
from typing import List, Tuple, Optional
import warnings
import tempfile
import os
import subprocess
import sys
import shutil
import gc
warnings.filterwarnings('ignore')

try:
    from ase.io import read, write
    from ase import Atoms
    from ase.data import covalent_radii
    from ase.data.colors import jmol_colors
    from ase.build import make_supercell
except ImportError as e:
    raise ImportError(f"ASE is required. Install with: pip install ase\nError: {e}")

# Ensure Qt platform plugins are discoverable for pip-installed PySide6/OVITO.
_pyside6_plugins = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir,
)
for _candidate in [
    os.path.join(sys.prefix, "lib", "python" + f"{sys.version_info.major}.{sys.version_info.minor}",
                 "site-packages", "PySide6", "Qt", "plugins"),
    os.path.join(sys.prefix, "lib", "qt6", "plugins"),
]:
    if os.path.isdir(_candidate):
        os.environ.setdefault("QT_PLUGIN_PATH", _candidate)
        break

# Optional imports for OVITO rendering (lazy-loaded to avoid startup instability).
OVITO_AVAILABLE = False
TACHYON_AVAILABLE = False
OSPRAY_AVAILABLE = False
OPENGL_AVAILABLE = False
ANARI_AVAILABLE = False
_OVITO_PROBED = False


def _probe_ovito_runtime() -> None:
    """Probe OVITO runtime capabilities lazily and cache the result."""
    global _OVITO_PROBED
    global OVITO_AVAILABLE, TACHYON_AVAILABLE, OSPRAY_AVAILABLE, OPENGL_AVAILABLE, ANARI_AVAILABLE
    if _OVITO_PROBED:
        return

    try:
        from ovito.io import import_file as _ovito_import_file
        from ovito.vis import Viewport as _ovito_viewport
        globals()["import_file"] = _ovito_import_file
        globals()["Viewport"] = _ovito_viewport
        OVITO_AVAILABLE = True
    except Exception:
        OVITO_AVAILABLE = False
        _OVITO_PROBED = True
        return

    try:
        from ovito.vis import TachyonRenderer as _tachyon_renderer
        globals()["TachyonRenderer"] = _tachyon_renderer
        TACHYON_AVAILABLE = True
    except Exception:
        TACHYON_AVAILABLE = False

    try:
        from ovito.vis import OSPRayRenderer as _ospray_renderer
        globals()["OSPRayRenderer"] = _ospray_renderer
        OSPRAY_AVAILABLE = True
    except Exception:
        OSPRAY_AVAILABLE = False

    try:
        from ovito.vis import OpenGLRenderer as _opengl_renderer
        globals()["OpenGLRenderer"] = _opengl_renderer
        OPENGL_AVAILABLE = True
    except Exception:
        OPENGL_AVAILABLE = False

    try:
        from ovito.vis import AnariRenderer as _anari_renderer
        globals()["AnariRenderer"] = _anari_renderer
        ANARI_AVAILABLE = True
    except Exception:
        ANARI_AVAILABLE = False

    _OVITO_PROBED = True


class CatalystSurfaceVisualizer:
    """High-quality visualizer for catalyst surfaces using OVITO renderers."""

    def __init__(
        self,
        elevation: float = 30,
        azimuth: float = 45,
        scale: float = 1.0,
        quality: str = 'high',
        background_color: str = 'white',
        color_scheme: str = 'jmol',
        renderer: str = 'tachyon',
        show_cell: bool = False,
        auto_expand: bool = True,
        expand_threshold: int = 15,
        supercell: Tuple[int, int, int] = (2, 2, 1)
    ):
        """
        Initialize the visualizer.

        Args:
            elevation: Viewing angle elevation in degrees (default: 30)
            azimuth: Viewing angle azimuth in degrees (default: 45)
            scale: Scale factor for atom sizes (default: 1.0)
            quality: Rendering quality ('low', 'medium', 'high', 'ultra')
            background_color: Background color
            color_scheme: Color scheme ('jmol', 'cpk')
            renderer: Rendering backend ('tachyon', 'ospray', 'opengl', 'anari')
            show_cell: Show unit cell boundaries (default: False)
            auto_expand: Automatically expand small structures (default: True)
            expand_threshold: Minimum number of atoms before expansion (default: 15)
            supercell: Supercell dimensions for expansion (default: 2x2x1)
        """
        self.elevation = elevation
        self.azimuth = azimuth
        self.scale = scale
        self.quality = quality
        self.background_color = background_color
        self.color_scheme = color_scheme
        self.renderer = renderer
        self.show_cell = show_cell
        self.auto_expand = auto_expand
        self.expand_threshold = expand_threshold
        self.supercell = supercell

        # Check renderer availability
        ovito_renderers = {
            'tachyon': TACHYON_AVAILABLE,
            'ospray': OSPRAY_AVAILABLE,
            'opengl': OPENGL_AVAILABLE,
            'anari': ANARI_AVAILABLE
        }

        if renderer in ovito_renderers:
            use_direct = os.environ.get("CATDT_OVITO_DIRECT", "0") == "1"
            if use_direct:
                _probe_ovito_runtime()
                if not OVITO_AVAILABLE:
                    raise RuntimeError(
                        "OVITO Python module is required for direct rendering but is unavailable."
                    )
                ovito_renderers = {
                    'tachyon': TACHYON_AVAILABLE,
                    'ospray': OSPRAY_AVAILABLE,
                    'opengl': OPENGL_AVAILABLE,
                    'anari': ANARI_AVAILABLE
                }
                if not ovito_renderers.get(renderer, False):
                    if TACHYON_AVAILABLE:
                        warnings.warn(f"{renderer.upper()} renderer unavailable, fallback to TACHYON.")
                        self.renderer = 'tachyon'
                    else:
                        raise RuntimeError(
                            f"Renderer {renderer} unavailable and TACHYON fallback is not available in this OVITO build"
                        )

        # Quality settings
        quality_settings = {
            'low': {
                'dpi': 150, 'resolution': 40, 'strides': 2,
                'ovito_width': 800, 'ovito_height': 800,
                'samples_per_pixel': 4, 'aa_level': 2
            },
            'medium': {
                'dpi': 200, 'resolution': 60, 'strides': 1,
                'ovito_width': 1200, 'ovito_height': 1200,
                'samples_per_pixel': 8, 'aa_level': 4
            },
            'high': {
                'dpi': 300, 'resolution': 100, 'strides': 1,
                'ovito_width': 2048, 'ovito_height': 2048,
                'samples_per_pixel': 16, 'aa_level': 8
            },
            'ultra': {
                'dpi': 600, 'resolution': 150, 'strides': 1,
                'ovito_width': 4096, 'ovito_height': 4096,
                'samples_per_pixel': 32, 'aa_level': 16
            }
        }

        settings = quality_settings.get(quality, quality_settings['high'])
        self.dpi = settings['dpi']
        self.resolution = settings['resolution']
        self.strides = settings['strides']
        self.ovito_width = settings['ovito_width']
        self.ovito_height = settings['ovito_height']
        self.samples_per_pixel = settings['samples_per_pixel']
        self.aa_level = settings['aa_level']

    def expand_structure_if_needed(self, atoms: Atoms) -> Tuple[Atoms, bool]:
        """
        Expand structure if it has fewer atoms than threshold.

        Args:
            atoms: ASE Atoms object

        Returns:
            (expanded_atoms, was_expanded): Tuple of atoms and whether expansion occurred
        """
        if not self.auto_expand:
            return atoms, False

        if len(atoms) < self.expand_threshold:
            print(f"Structure has {len(atoms)} atoms (< {self.expand_threshold})")
            print(f"Expanding to {self.supercell[0]}x{self.supercell[1]}x{self.supercell[2]} supercell...")

            # Create supercell matrix
            P = np.diag(self.supercell)
            expanded = make_supercell(atoms, P)

            print(f"Expanded structure now has {len(expanded)} atoms")
            return expanded, True

        return atoms, False

    def get_atom_color(self, atomic_number: int) -> np.ndarray:
        """Get color for atomic number from jmol scheme."""
        return np.array(jmol_colors[atomic_number])

    def get_atom_radius(self, atomic_number: int) -> float:
        """Get radius for atomic number."""
        return covalent_radii[atomic_number] * self.scale

    def render_with_ovito(
        self,
        atoms: Atoms,
        output_file: str,
        renderer_type: str = 'tachyon'
    ) -> None:
        """
        Render structure using OVITO renderers.

        By default, rendering is isolated in a subprocess to avoid hard
        process crashes (segfault) from OVITO/Tachyon C++ backends.
        """
        use_direct = os.environ.get("CATDT_OVITO_DIRECT", "0") == "1"
        if use_direct:
            self._render_with_ovito_direct(atoms, output_file, renderer_type)
        else:
            self._render_with_ovito_subprocess(atoms, output_file, renderer_type)

    def _render_with_ovito_subprocess(
        self,
        atoms: Atoms,
        output_file: str,
        renderer_type: str = 'tachyon'
    ) -> None:
        """Isolate OVITO rendering in a subprocess for robustness."""
        temp_fd, temp_file = tempfile.mkstemp(suffix='.xyz')
        os.close(temp_fd)

        try:
            write(temp_file, atoms)

            ovitos_exec = os.environ.get("CATDT_OVITOS_EXECUTABLE", "").strip() or shutil.which("ovitos")
            interpreter = ovitos_exec if ovitos_exec else sys.executable
            cmd = [
                interpreter,
                str(Path(__file__).resolve()),
                temp_file,
                "-o",
                output_file,
                "--renderer",
                renderer_type,
                "--quality",
                self.quality,
                "--elevation",
                str(self.elevation),
                "--azimuth",
                str(self.azimuth),
                "--scale",
                str(self.scale),
                "--background",
                str(self.background_color),
                "--color-scheme",
                str(self.color_scheme),
                "--expand-threshold",
                str(self.expand_threshold),
                "--supercell",
                f"{self.supercell[0]},{self.supercell[1]},{self.supercell[2]}",
            ]
            if self.show_cell:
                cmd.append("--show-cell")
            if not self.auto_expand:
                cmd.append("--no-expand")

            env = os.environ.copy()
            env["CATDT_OVITO_DIRECT"] = "1"
            timeout_s = int(str(os.environ.get("CATDT_OVITO_RENDER_TIMEOUT_SEC", "180")).strip() or "180")
            try:
                result = subprocess.run(
                    cmd,
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=max(10, timeout_s),
                )
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(
                    f"OVITO subprocess render timeout after {max(10, timeout_s)}s. "
                    f"cmd={' '.join(cmd)}"
                ) from exc
            if result.returncode != 0:
                err = (result.stderr or result.stdout or "").strip()
                if len(err) > 2000:
                    err = err[-2000:]
                raise RuntimeError(
                    f"OVITO subprocess render failed (code={result.returncode}): {err}"
                )
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def _render_with_ovito_direct(
        self,
        atoms: Atoms,
        output_file: str,
        renderer_type: str = 'tachyon'
    ) -> None:
        """
        Render structure using OVITO renderers.

        Args:
            atoms: ASE Atoms object
            output_file: Output filename
            renderer_type: 'tachyon', 'ospray', 'opengl', or 'anari'
        """
        from ovito import scene as ovito_scene
        from ovito.io import import_file
        from ovito.vis import Viewport

        # Save atoms to temporary file for OVITO to import
        temp_fd, temp_file = tempfile.mkstemp(suffix='.xyz')
        os.close(temp_fd)
        pipeline = None
        data = None
        viewport = None
        renderer = None

        try:
            # Write atoms to temporary xyz file
            write(temp_file, atoms)

            # Import with OVITO
            pipeline = import_file(temp_file)

            # Add pipeline to scene for rendering
            pipeline.add_to_scene()

            # Compute pipeline
            data = pipeline.compute()

            # Configure cell display
            if data.cell:
                try:
                    # Control cell visualization based on user preference
                    data.cell.vis.enabled = self.show_cell
                except:
                    pass  # Cell visualization may not be available

            # Create viewport
            viewport = Viewport()
            viewport.type = Viewport.Type.Perspective

            # Set camera position
            positions = atoms.get_positions()
            center = positions.mean(axis=0)
            bbox_size = positions.max(axis=0) - positions.min(axis=0)
            max_dim = np.max(bbox_size)
            camera_distance = max_dim * 3.0

            # Camera position in spherical coordinates (30° elevation, looking down)
            elev_rad = np.radians(self.elevation)
            azim_rad = np.radians(self.azimuth)
            cam_x = center[0] + camera_distance * np.cos(elev_rad) * np.cos(azim_rad)
            cam_y = center[1] + camera_distance * np.cos(elev_rad) * np.sin(azim_rad)
            cam_z = center[2] + camera_distance * np.sin(elev_rad)

            viewport.camera_pos = (cam_x, cam_y, cam_z)
            viewport.camera_dir = (center[0] - cam_x, center[1] - cam_y, center[2] - cam_z)
            viewport.fov = np.radians(35)

            # Create renderer based on type
            if renderer_type == 'tachyon':
                from ovito.vis import TachyonRenderer
                renderer = TachyonRenderer()
                renderer.antialiasing = True
                renderer.direct_light = True
                renderer.shadows = True
                renderer.ambient_occlusion = True
                renderer.ambient_occlusion_samples = self.samples_per_pixel
                renderer.direct_light_intensity = 0.9

            elif renderer_type == 'ospray':
                from ovito.vis import OSPRayRenderer
                renderer = OSPRayRenderer()
                renderer.refinement_iterations = self.samples_per_pixel
                renderer.denoising_enabled = True
                renderer.ambient_light_enabled = True
                renderer.direct_light_enabled = True
                renderer.direct_light_intensity = 1.0
                renderer.ambient_brightness = 0.5

            elif renderer_type == 'opengl':
                from ovito.vis import OpenGLRenderer
                renderer = OpenGLRenderer()
                renderer.antialiasing_level = self.aa_level

            elif renderer_type == 'anari':
                from ovito.vis import AnariRenderer
                renderer = AnariRenderer()
                try:
                    renderer.denoising_enabled = True
                    renderer.ambient_occlusion_samples = self.samples_per_pixel
                except AttributeError:
                    pass

            else:
                raise ValueError(f"Unknown renderer type: {renderer_type}")

            # Render to file
            viewport.render_image(
                filename=output_file,
                size=(self.ovito_width, self.ovito_height),
                renderer=renderer,
                crop=False,
                alpha=False
            )

        finally:
            # Clean up temporary file
            if os.path.exists(temp_file):
                os.unlink(temp_file)

            # Remove from scene
            try:
                if pipeline is not None:
                    pipeline.remove_from_scene()
            except Exception:
                pass
            try:
                ovito_scene.pipelines.clear()
            except Exception:
                pass
            del renderer
            del viewport
            del data
            del pipeline
            gc.collect()

    def visualize_structure(
        self,
        atoms: Atoms,
        output_file: Optional[str] = None,
        title: Optional[str] = None,
        show: bool = False
    ) -> None:
        """Visualize a single structure."""
        _ = title
        _ = show
        # Expand if needed
        atoms, was_expanded = self.expand_structure_if_needed(atoms)

        if output_file is None:
            output_file = 'output.png'
        print(f"Rendering with {self.renderer.upper()}...")
        self.render_with_ovito(atoms, output_file, self.renderer)
        print(f"✓ Saved: {output_file}")

    def visualize_trajectory(
        self,
        trajectory: List[Atoms],
        output_file: str,
        fps: int = 2,
        titles: Optional[List[str]] = None,
        show_progress: bool = True
    ) -> None:
        """Visualize trajectory as GIF animation."""
        if not output_file.endswith('.gif'):
            output_file = str(Path(output_file).with_suffix('.gif'))

        n_frames = len(trajectory)
        if titles is None:
            titles = [f"Frame {i+1}/{n_frames}" for i in range(n_frames)]

        # Check if we need to expand structures
        if self.auto_expand and len(trajectory[0]) < self.expand_threshold:
            print(f"Trajectory frames have {len(trajectory[0])} atoms (< {self.expand_threshold})")
            print(f"Expanding all frames to {self.supercell[0]}x{self.supercell[1]}x{self.supercell[2]} supercell...")
            P = np.diag(self.supercell)
            trajectory = [make_supercell(atoms, P) for atoms in trajectory]
            print(f"Expanded frames now have {len(trajectory[0])} atoms")

        print(f"Rendering trajectory with {self.renderer.upper()}...")
        temp_dir = tempfile.mkdtemp()
        frame_files = []

        for i, atoms in enumerate(trajectory):
            if show_progress:
                print(f"\rRendering frame {i + 1}/{n_frames}...", end='', flush=True)
            frame_file = os.path.join(temp_dir, f"frame_{i:04d}.png")
            self.render_with_ovito(atoms, frame_file, self.renderer)
            frame_files.append(frame_file)

        from PIL import Image
        images = [Image.open(f) for f in frame_files]
        try:
            images[0].save(
                output_file,
                save_all=True,
                append_images=images[1:],
                duration=int(1000 / fps),
                loop=0,
            )
        finally:
            for image in images:
                try:
                    image.close()
                except Exception:
                    pass

        import shutil
        shutil.rmtree(temp_dir)
        gc.collect()

        if show_progress:
            print(f"\n✓ Saved: {output_file}")
            print(f"  Frames: {n_frames}, FPS: {fps}, Duration: {n_frames/fps:.1f}s")


def auto_detect_file_type(filename: str) -> Tuple[bool, int]:
    """Auto-detect if file is trajectory or single structure."""
    try:
        structures = read(filename, index=':')
        if isinstance(structures, list):
            return len(structures) > 1, len(structures)
        else:
            return False, 1
    except:
        try:
            read(filename)
            return False, 1
        except:
            raise ValueError(f"Cannot read file: {filename}")


def main():
    """Main CLI."""
    parser = argparse.ArgumentParser(
        description='High-quality catalyst surface visualizer (OVITO renderers only)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single structure with Tachyon (default)
  %(prog)s structure.xyz

  # Render with OSPRay ray tracer
  %(prog)s structure.xyz --renderer ospray

  # Ultra quality with OSPRay
  %(prog)s structure.xyz --renderer ospray --quality ultra

  # Show unit cell boundaries
  %(prog)s structure.xyz --show-cell

  # Disable auto-expansion for small structures
  %(prog)s structure.xyz --no-expand

  # Trajectory animation
  %(prog)s trajectory.traj --fps 2

  # Custom angle and larger atoms
  %(prog)s structure.xyz --elevation 60 --azimuth 30 --scale 1.5

Rendering Backends:
  tachyon:    Ray tracer with high-quality shadows [DEFAULT]
  ospray:     Intel's ray tracer with advanced rendering
  opengl:     Fast OpenGL renderer for real-time rendering
  anari:      Hardware accelerated ray-tracing (NVIDIA VisRTX)
Installation:
  pip install ase pillow
  pip install ovito  # For tachyon, ospray, opengl, and anari
        """
    )

    parser.add_argument('input', help='Input file (xyz, traj, cif, POSCAR, etc.)')
    parser.add_argument('-o', '--output', help='Output file', default=None)
    parser.add_argument('--renderer', choices=['tachyon', 'ospray', 'opengl', 'anari'],
                       default='tachyon', help='Rendering backend (default: tachyon)')
    parser.add_argument('--elevation', type=float, default=30, help='View elevation (default: 30°)')
    parser.add_argument('--azimuth', type=float, default=45, help='View azimuth (default: 45°)')
    parser.add_argument('--scale', type=float, default=1.0, help='Atom size scale (default: 1.0)')
    parser.add_argument('--quality', choices=['low', 'medium', 'high', 'ultra'],
                       default='high', help='Rendering quality (default: high)')
    parser.add_argument('--fps', type=int, default=2, help='GIF frames per second (default: 2)')
    parser.add_argument('--background', default='white', help='Background color (default: white)')
    parser.add_argument('--color-scheme', choices=['jmol', 'cpk'], default='jmol',
                       help='Color scheme (default: jmol)')
    parser.add_argument('--show-cell', action='store_true', help='Show unit cell boundaries')
    parser.add_argument('--no-expand', action='store_true', help='Disable auto-expansion for small structures')
    parser.add_argument('--expand-threshold', type=int, default=15,
                       help='Minimum atoms before expansion (default: 15)')
    parser.add_argument('--supercell', type=str, default='2,2,1',
                       help='Supercell dimensions as "nx,ny,nz" (default: 2,2,1)')
    parser.add_argument('--show', action='store_true', help='Display the visualization')

    args = parser.parse_args()

    # Check input
    if not Path(args.input).exists():
        print(f"Error: File not found: {args.input}")
        return 1

    # Check OVITO availability
    if args.renderer in ['tachyon', 'ospray', 'opengl', 'anari']:
        _probe_ovito_runtime()
        ovitos_exec = os.environ.get("CATDT_OVITOS_EXECUTABLE", "").strip() or shutil.which("ovitos")
        if not OVITO_AVAILABLE and not ovitos_exec:
            print(f"Error: {args.renderer} renderer requires OVITO runtime.")
            print("Install OVITO Python module or provide `ovitos` executable in PATH.")
            return 1

    # Parse supercell dimensions
    try:
        supercell = tuple(map(int, args.supercell.split(',')))
        if len(supercell) != 3:
            raise ValueError
    except:
        print("Error: --supercell must be in format 'nx,ny,nz' (e.g., '2,2,1')")
        return 1

    # Detect file type
    is_trajectory, n_frames = auto_detect_file_type(args.input)

    # Create visualizer
    visualizer = CatalystSurfaceVisualizer(
        elevation=args.elevation,
        azimuth=args.azimuth,
        scale=args.scale,
        quality=args.quality,
        background_color=args.background,
        color_scheme=args.color_scheme,
        renderer=args.renderer,
        show_cell=args.show_cell,
        auto_expand=not args.no_expand,
        expand_threshold=args.expand_threshold,
        supercell=supercell
    )

    # Output filename
    if args.output is None:
        input_path = Path(args.input)
        args.output = input_path.with_suffix('.gif' if is_trajectory else '.png')

    # Info
    print("=" * 70)
    print("Catalyst Surface Visualizer")
    print("=" * 70)
    print(f"Input: {args.input}")
    print(f"Renderer: {args.renderer}")
    print(f"Type: {'Trajectory' if is_trajectory else 'Single structure'}")
    if is_trajectory:
        print(f"Frames: {n_frames}")
    print(f"Quality: {args.quality}")
    if args.renderer in ['tachyon', 'ospray', 'opengl', 'anari']:
        print(f"Resolution: {visualizer.ovito_width}x{visualizer.ovito_height}")
        if args.renderer != 'opengl':
            print(f"Samples per pixel: {visualizer.samples_per_pixel}")
    print(f"View: elevation={args.elevation}°, azimuth={args.azimuth}°")
    print(f"Show cell: {args.show_cell}")
    print(f"Auto-expand: {not args.no_expand} (threshold: {args.expand_threshold} atoms)")
    if not args.no_expand:
        print(f"Supercell: {supercell[0]}x{supercell[1]}x{supercell[2]}")
    print("=" * 70)

    if is_trajectory:
        trajectory = read(args.input, index=':')
        titles = [f"Frame {i+1}/{n_frames}" for i in range(n_frames)]
        visualizer.visualize_trajectory(trajectory, str(args.output), args.fps, titles)
    else:
        atoms = read(args.input)
        print(f"Atoms: {len(atoms)}")
        print(f"Elements: {', '.join(sorted(set(atoms.get_chemical_symbols())))}")
        visualizer.visualize_structure(atoms, str(args.output), show=args.show)

    print("=" * 70)
    print("Done!")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    exit(main())
