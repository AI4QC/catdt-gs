"""Performs sampling of surface reconstructions using an MCMC-based algorithm."""

import logging
import random
from collections import defaultdict
from collections.abc import Iterable
from collections import Counter
from pathlib import Path

import numpy as np

from mcmc.events.criterion import (
    DistanceCriterion,
    MetropolisCriterion,
    TestingCriterion,
)
from mcmc.events.event import Change, Exchange
from mcmc.slab import change_site
from mcmc.events.proposal import ChangeProposal, SwitchProposal
from mcmc.system import SurfaceSystem
from mcmc.utils import create_anneal_schedule, setup_folders
from mcmc.utils.clustering import (
    find_closest_points_indices,
    get_cluster_centers,
)
from mcmc.utils.plot import plot_clustering_results


class MCMC:
    """MCMC-based class for sampling surface reconstructions."""

    def __init__(
        self,
        adsorbates=None,
        canonical=False,
        num_ads_atoms=0,
        adsorbate_counts=None,
        testing=False,
        filter_distance: float = 0.0,
        **kwargs,
    ) -> None:
        """Initialize the MCMC class.

        Args:
            adsorbates (list, optional): The list of adsorbates. Defaults to None.
            canonical (bool, optional): If True, perform canonical sampling. Defaults to False.
            num_ads_atoms (int, optional): The number of adsorbed atoms. Defaults to 0.
            testing (bool, optional): If True, perform testing. Defaults to False.
            filter_distance (float, optional): The distance for filtering. Defaults to 0.0.
            **kwargs: Additional keyword arguments.

        Attributes:
            canonical (bool): If True, perform canonical sampling.
            num_ads_atoms (int): The number of adsorbed atoms.
            testing (bool): If True, perform testing.
            adsorbates (list): The list of adsorbates.
            filter_distance (float): The distance for filtering.
            kwargs: Additional keyword arguments.
            surface: The surface system. Defaults to None.
            total_sweeps (int): The total number of sweeps. Defaults to 100.
            sweep_size (int): The number of steps to perform in each sweep. Defaults to 100.
            temp (float): The temperature parameter used in the Metropolis-Hastings algorithm for MC
                simulations. Defaults to 1.0.
            alpha (float): The alpha parameter used in the annealing schedule. Defaults to 1.0.
            logger (logging.Logger): The logger object. Defaults to None.
            run_folder (str): The folder in which to save the results of the simulation. Defaults to
                "".

        Raises:
            AssertionError: If canonical sampling is selected but the number of adsorbed atoms is
            fewer than 0.
        """
        self.adsorbates = adsorbates
        self.canonical = canonical
        self.num_ads_atoms = num_ads_atoms
        self.adsorbate_counts = dict(adsorbate_counts or {})
        self.testing = testing
        self.filter_distance = filter_distance
        self.kwargs = kwargs

        # Initialize here for subsequent runs
        self.surface = None
        self.total_sweeps = 100
        self.sweep_size = 100
        self.temp = 1.0
        self.alpha = 1.0
        self.logger = None
        self.run_folder = ""

        if self.canonical:
            if self.adsorbate_counts:
                self.num_ads_atoms = sum(int(v) for v in self.adsorbate_counts.values())
            assert self.num_ads_atoms > 0, (
                "for canonical runs, need number of adsorbed atoms greater than 0"
            )

    def initialize(
        self,
        even_adsorption_sites: bool = False,
        perform_annealing: bool = False,
        anneal_schedule: list | np.ndarray | None = None,
        multiple_anneal: bool = False,
        run_folder: Path | str | None = None,
        **kwargs,
    ) -> np.ndarray:
        """Initialize the MCMC simulation by setting up the run folder, preparing the canonical
        slabs, and creating the annealing schedule.

        Args:
            even_adsorption_sites (bool, optional): If True, evenly adsorb the sites. Defaults to
                False.
            perform_annealing (bool, optional): If True, perform annealing. Defaults to False.
            anneal_schedule (list | np.ndarray, optional): The annealing schedule. Defaults to None.
            multiple_anneal (bool, optional): If True, perform multiple annealing. Defaults to
                False.
            run_folder (Path | str | None, optional): The folder in which to save the results of the
                simulation. Defaults to None.
            **kwargs: Additional keyword arguments.

        Returns:
            np.ndarray: The annealing schedule.
        """
        if not run_folder:
            self.run_folder = setup_folders(
                self.surface.surface_name,
                canonical=self.canonical,
                total_sweeps=self.total_sweeps,
                start_temp=self.temp,
                alpha=self.alpha,
                **kwargs,
            )
        else:
            self.run_folder = Path(run_folder)

        self.logger.info("Using run folder %s", self.run_folder)
        self.run_folder.mkdir(parents=True, exist_ok=True)

        if self.canonical:
            self.prepare_canonical(even_adsorption_sites=even_adsorption_sites)

        if isinstance(anneal_schedule, Iterable):
            temp_list = anneal_schedule  # user-defined annealing schedule
        elif perform_annealing:
            temp_list = create_anneal_schedule(
                start_temp=self.temp,
                total_sweeps=self.total_sweeps,
                alpha=self.alpha,
                multiple_anneal=multiple_anneal,
                save_folder=self.run_folder,
            )
        else:
            temp_list = np.repeat(self.temp, self.total_sweeps)  # constant temperature
        return temp_list

    def prepare_canonical(self, even_adsorption_sites: bool = False):
        """Prepare a canonical slab by performing semi-grand canonical adsorption runs until the
        desired number of adsorbed atoms are obtained.

        Args:
            even_adsorption_sites (bool, optional): If True, evenly adsorb the sites. Defaults to
            False.

        Raises:
            AssertionError: If the number of adsorbed atoms is less than 0.
        """
        assert (
            self.num_ads_atoms > 0
        ), "for canonical runs, need number of adsorbed atoms greater than 0"

        if even_adsorption_sites:
            self.logger.info("Evenly adsorbing sites")
            # Do clustering
            centers, labels = get_cluster_centers(
                self.surface.ads_coords[:, :2], self.num_ads_atoms
            )
            sites_idx = find_closest_points_indices(self.surface.ads_coords[:, :2], centers, labels)
            plot_clustering_results(
                self.surface.ads_coords,
                self.num_ads_atoms,
                labels,
                sites_idx,
                save_folder=self.run_folder,
            )

            for site_idx in sites_idx:
                self.step_semigrand(site_idx=site_idx)
        elif self.adsorbate_counts:
            self.logger.info("Randomly adsorbing fixed composition: %s", self.adsorbate_counts)
            n_sites = len(self.surface.ads_coords)
            target_counts = {
                str(elem): int(count)
                for elem, count in self.adsorbate_counts.items()
                if int(count) > 0
            }
            target_total = sum(target_counts.values())
            if target_total > n_sites:
                raise ValueError(
                    f"adsorbate_counts target {target_total} exceeds n_sites={n_sites}: "
                    f"{target_counts}"
                )
            site_indices = list(range(n_sites))
            random.shuffle(site_indices)
            cursor = 0
            for elem, count in target_counts.items():
                if elem not in self.adsorbates:
                    raise ValueError(
                        f"adsorbate_counts contains {elem!r}, but adsorbates={self.adsorbates}"
                    )
                for site_idx in site_indices[cursor : cursor + count]:
                    self.surface = change_site(self.surface, site_idx, elem)
                cursor += count
            self.num_ads_atoms = target_total
        else:
            self.logger.info("Randomly adsorbing sites")
            # Cap the target at the number of available adsorption sites —
            # otherwise, if num_ads_atoms > len(ads_coords), the loop can
            # never reach the target and runs forever (each site holds at
            # most one adsorbate). Seen on post-MD Cu(111) slabs where
            # AdsorbateSiteFinder yields few sites.
            n_sites = len(self.surface.ads_coords)
            if self.num_ads_atoms > n_sites:
                self.logger.warning(
                    "num_ads_atoms=%s exceeds n_sites=%s; capping target",
                    self.num_ads_atoms, n_sites,
                )
                self.num_ads_atoms = n_sites
            # Safety bound: prepare_canonical shouldn't need more than
            # 20 * n_sites rejections to reach target on a healthy system.
            max_attempts = 20 * max(n_sites, 1)
            attempts = 0
            while self.surface.num_adsorbates < self.num_ads_atoms:
                self.step_semigrand()
                attempts += 1
                if attempts > max_attempts:
                    self.logger.warning(
                        "prepare_canonical: giving up after %d attempts "
                        "(reached %d/%d adsorbates)",
                        attempts, self.surface.num_adsorbates, self.num_ads_atoms,
                    )
                    break

        self.surface.real_atoms.write(
            self.run_folder / f"{self.surface.surface_name}_canonical_init.cif"
        )

    def step_canonical(self, iter_num: int = 1) -> bool:
        """Performs canonical sampling. Switches the adsorption sites of two
        adsorbates and calculates acceptance based on specified criterion.

        Args:
            iter_num (int, optional): The iteration number of the simulation. Defaults to 1.

        Returns:
            bool: Whether the proposed change was accepted or not.
        """
        if iter_num % self.sweep_size == 0:
            self.logger.info("At iter %s", iter_num)
            plot_specific_distance_weights = True
        else:
            plot_specific_distance_weights = False

        proposal = SwitchProposal(
            system=self.surface,
            adsorbate_list=self.adsorbates.copy(),
            require_per_atom_energies=self.kwargs.get("require_per_atom_energies", False),
            require_distance_decay=self.kwargs.get("require_distance_decay", False),
            temp=self.temp,
            run_folder=self.run_folder,
            plot_specific_distance_weights=plot_specific_distance_weights,
            run_iter=iter_num,
        )
        self.logger.debug("\n we are at iter %s", iter_num)

        if self.filter_distance > 0:
            criterion = DistanceCriterion(filter_distance=self.filter_distance)
            self.logger.debug("Using distance filter")
        elif self.testing:
            criterion = TestingCriterion()
            self.logger.debug("Using test criterion, always accept")
        else:
            criterion = MetropolisCriterion(self.temp)
            self.logger.debug("Using Metropolis criterion")

        event = Exchange(self.surface, proposal, criterion)

        accept, self.surface = event.acceptance()
        return accept

    def step_semigrand(self, iter_num: int = 1, site_idx: int | None = None) -> bool:
        """Performs a semigrand canonical sampling iteration. It randomly chooses a site to change
        identity in a slab, adds or removes an atom from the site, optionally performs relaxation,
        calculates acceptance based on specified criterion.

        Args:
            iter_num (int, optional): The iteration number of the simulation.
            site_idx (int, optional): Specify the index of the site to switch.

        Returns:
            bool: Whether the proposed change was accepted or not.
        """
        self.logger.debug("\n we are at iter %s", iter_num)

        proposal = ChangeProposal(
            system=self.surface,
            adsorbate_list=self.adsorbates.copy(),
            site_idx=site_idx,
        )

        if self.filter_distance > 0:
            criterion = DistanceCriterion(filter_distance=self.filter_distance)
            self.logger.debug("Using distance filter")
        elif self.testing:
            criterion = TestingCriterion()
            self.logger.debug("Using test criterion, always accept")
        else:
            criterion = MetropolisCriterion(self.temp)
            self.logger.debug("Using Metropolis criterion")

        event = Change(self.surface, proposal, criterion)

        accept, self.surface = event.acceptance()
        return accept

    def sweep(self, i: int = 0) -> dict:
        """Perform MC sweep.

        Args:
            i (int, optional): The sweep number. Defaults to 0.

        Returns:
            dict: A dictionary containing the history, trajectory, energy, adsorption count, and
                acceptance rate.
        """
        num_accept = 0
        self.logger.info("In sweep %s out of %s", i + 1, self.total_sweeps)
        for j in range(self.sweep_size):
            run_idx = self.sweep_size * i + j + 1
            if self.canonical:
                accept = self.step_canonical(iter_num=run_idx)
            else:
                accept = self.step_semigrand(iter_num=run_idx)
            num_accept += accept

        # Save structure and traj for easy viewing
        self.surface.save_structures(sweep_num=i + 1, save_folder=self.run_folder)

        surface = self.surface.copy(copy_calc=False)
        surface.unset_calc()
        return {
            "history": surface,
            "trajectory": self.surface.relax_traj,
            "energy": self.surface.get_surface_energy(),
            "adsorption_count": self.surface.num_adsorbates,
            "acceptance_rate": num_accept / self.sweep_size,
        }

    def run(
        self,
        surface: SurfaceSystem,
        logger: logging.Logger | None = None,
        total_sweeps: int = 100,
        sweep_size: int = 20,
        start_temp: float = 1.0,
        perform_annealing: bool = True,
        alpha: float = 0.99,
        multiple_anneal: bool = False,
        anneal_schedule: list | None = None,
        run_folder: str | None = None,
        starting_iteration: list = 0,
        even_adsorption_sites: bool = False,
        **kwargs,
    ) -> dict:
        """This function runs an MC simulation for a given number of sweeps and temperature, and
        returns the history of the simulation along with summary statistics.

        Args:
            surface (SurfaceSystem): The surface system on which the MCMC simulation is to be
                performed.
            logger (logging.Logger, optional): The logger object. Defaults to None.
            total_sweeps (int, optional): The number of MCMC sweeps to perform. Defaults to 100.
            sweep_size : int, optional
                The number of steps to perform in each sweep. Defaults to 20.
            start_temp : float, optional
                The temperature parameter is used in the Metropolis-Hastings algorithm for MC
                simulations. It controls the probability of accepting a proposed move during the
                simulation. A higher temperature leads to a higher probability of accepting a move,
                while a lower temperature leads to a lower probability of accepting a move. Defaults
                to 1.0.
            perform_annealing : bool, optional
                If True, perform annealing. Defaults to True.
            alpha : float, optional
                The alpha parameter is a value between 0 and 1 that determines the annealing rate.
                A higher alpha results in a slower annealing rate, while a lower alpha results in
                a faster annealing rate. Defaults to 0.99.
            multiple_anneal : bool, optional
                If True, perform multiple annealing. Defaults to False.
            anneal_schedule : list, optional
                The annealing schedule. Defaults to None.
            run_folder : str, optional
                The folder in which to save the results of the simulation. Defaults to None.
            starting_iteration : list, optional
                The starting iteration number of the simulation. Defaults to 0.
            even_adsorption_sites : bool, optional
                If True, evenly adsorb the sites. Defaults to False.
            **kwargs : dict, optional

        Returns:
            dict: A dictionary containing the history, trajectory, energy, adsorption count, and
                acceptance rate.
        """
        self.surface = surface
        self.logger = logger or logging.getLogger(__name__)
        self.total_sweeps = total_sweeps
        self.sweep_size = sweep_size
        self.temp = start_temp
        self.alpha = alpha

        temp_list = self.initialize(
            even_adsorption_sites=even_adsorption_sites,
            perform_annealing=perform_annealing,
            anneal_schedule=anneal_schedule,
            multiple_anneal=multiple_anneal,
            run_folder=run_folder,
        )
        self.logger.info("There are %d atoms in pristine slab", self.surface.num_pristine_atoms)
        self.logger.info(
            "Running with total_sweeps = %d, sweep_size = %d, start_temp = %.3f",
            total_sweeps,
            sweep_size,
            start_temp,
        )
        self.logger.info("Starting with iteration %d", starting_iteration)
        self.logger.info("Temperature schedule is: %s", [f"{temp:.3f}" for temp in temp_list])

        # Perform MC sweeps
        results = defaultdict(list)  # TODO: make a dataclass
        for i in range(starting_iteration, self.total_sweeps):
            self.temp = temp_list[i]
            sweep_result = self.sweep(i=i)  # TODO change to .step
            results["history"].append(sweep_result["history"])
            results["trajectories"].append(sweep_result["trajectory"])
            results["energy_hist"].append(sweep_result["energy"])
            results["frac_accept_hist"].append(sweep_result["acceptance_rate"])
            results["adsorption_count_hist"].append(sweep_result["adsorption_count"])

        return results


if __name__ == "__main__":
    pass
