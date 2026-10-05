import os

from byteff2.train import FFTrainer
from bytemol.utils import setup_default_logging

# Showcase only: fine-tunes optimal.pt on example.h5 + the 60-geometry Cl- pilot set.
CONFIG = "train_cl_pilot.yaml"

os.chdir(os.path.dirname(os.path.abspath(__file__)))
logger = setup_default_logging()
trainer = FFTrainer(CONFIG, timestamp=False, restart=False)
trainer.train_loop()
logger.info("Training finished!")
