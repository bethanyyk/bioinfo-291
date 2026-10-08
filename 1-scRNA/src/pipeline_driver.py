"""
Stage-based pipeline driver with skip logic, --force and --stage flags.

Orchestrates a sequence of analysis stages, each of which:
  - Has an input/output contract (what data it expects, what it produces)
  - Can be run in isolation (e.g., re-running abundance without re-clustering)
  - Automatically skips if its outputs exist (unless --force is passed)
  - Can be selected individually (--stage abundance)

This decouples stage implementation from orchestration, making the pipeline
reproducible and debuggable.

Usage:
    driver = PipelineDriver(
        stages=[
            Stage('ingest', 'Load and parse raw data'),
            Stage('preprocess', 'QC and normalization'),
            Stage('annotate', 'Cell type assignment'),
        ],
        outputs={'ingest': ['data.h5ad'], 'preprocess': ['data.h5ad']},
    )
    driver.run(config, force=args.force, stage=args.stage)
"""

import logging
from pathlib import Path
from typing import Dict, List, Callable, Optional, Any
from dataclasses import dataclass


log = logging.getLogger(__name__)


@dataclass
class Stage:
    """A single analysis stage."""
    name: str
    description: str
    func: Optional[Callable] = None  # function to run; if None, skipped
    
    def __repr__(self):
        return f"{self.name}: {self.description}"


class PipelineDriver:
    """Orchestrates stages with skip logic and per-stage selection."""
    
    def __init__(
        self,
        stages: List[Stage],
        outputs: Dict[str, List[str]],
    ):
        """
        Initialize the pipeline.
        
        Parameters
        ----------
        stages : list of Stage
            Stages to run, in order.
        outputs : dict
            Maps stage name to list of output filenames (relative to results dir)
            that indicate the stage completed successfully.
        """
        self.stages = {s.name: s for s in stages}
        self.stage_order = [s.name for s in stages]
        self.outputs = outputs
    
    def _outputs_exist(self, stage_name: str, config: Any) -> bool:
        """Check if all outputs of a stage exist."""
        if stage_name not in self.outputs:
            return False
        
        results_dir = config.path('results')
        for output in self.outputs[stage_name]:
            p = results_dir / output
            if not p.exists():
                log.debug(f"  Output missing: {p.name}")
                return False
        return True
    
    def run(
        self,
        config: Any,
        force: bool = False,
        stage: Optional[str] = None,
        **stage_kwargs,
    ) -> Dict[str, bool]:
        """
        Run the pipeline.
        
        Parameters
        ----------
        config : Config
            Configuration object.
        force : bool
            If True, re-run all stages even if outputs exist.
        stage : str, optional
            If specified, run only this stage.
        **stage_kwargs
            Additional keyword arguments passed to each stage function.
        
        Returns
        -------
        dict
            Mapping of stage name to whether it was executed.
        """
        stages_to_run = [stage] if stage else self.stage_order
        results = {}
        
        for stage_name in stages_to_run:
            if stage_name not in self.stages:
                raise ValueError(f"Unknown stage: {stage_name}")
            
            s = self.stages[stage_name]
            outputs_exist = self._outputs_exist(stage_name, config)
            
            if outputs_exist and not force:
                log.info(f"[SKIP] {s.name}: outputs exist")
                results[stage_name] = False
                continue
            
            if s.func is None:
                log.info(f"[SKIP] {s.name}: no function defined")
                results[stage_name] = False
                continue
            
            log.info(f"[RUN]  {s.name}: {s.description}")
            try:
                s.func(config=config, **stage_kwargs)
                results[stage_name] = True
                log.info(f"[OK]   {s.name}")
            except Exception as e:
                log.error(f"[FAIL] {s.name}: {e}")
                raise
        
        return results
