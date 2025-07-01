# Project_Malaria

## Overview
Project_Malaria is a cheminformatics and machine learning platform designed for molecular generation, scoring, and optimization, with a focus on drug discovery for malaria and related research. The project leverages deep learning models, custom scoring functions, and flexible pipelines to accelerate molecular design and analysis.

## Features
- Deep generative models for molecule design (LibInvent, LinkInvent, Mol2Mol, etc.)
- Customizable scoring and filtering pipelines
- Support for transfer learning and staged learning
- Integration with RDKit and other cheminformatics tools
- Jupyter notebooks for interactive experimentation
- Modular plugin system for extending scoring and normalization

## Directory Structure
- `reinvent/` — Core source code for molecular generation and scoring
- `reinvent_plugins/` — Plugins for custom scoring, normalization, and components
- `configs/` — Example configuration files, scoring setups, and results
- `notebooks/` — Example and demo Jupyter notebooks
- `priors/` — Pretrained model weights and vocabularies
- `support/` — Utility scripts for data conversion and model management
- `tests/` — Unit and integration tests

## Getting Started
1. **Clone the repository:**
   ```sh
   git clone <repo-url>
   cd REINVENT4
   ```
2. **Install dependencies:**
   - Use `pip` and the provided `pyproject.toml` or `requirements.txt` (if available).
   - Set up a virtual environment for isolation.
3. **Run a demo notebook:**
   - Open a notebook from the `notebooks/` directory in Jupyter or VS Code.
   - Follow the instructions in the notebook to generate and score molecules.

## Usage
- Configure your experiment using TOML or JSON files in `configs/`.
- Run scripts or notebooks to train, sample, or score molecules.
- Extend functionality by adding plugins in `reinvent_plugins/`.

## Configuration Files
All configuration files for this project were written by Andrew Tom Mathew (BL.EN.U4CSE23269).

## License
This project is licensed under the terms of the LICENSE file.

This project uses or modifies code from the REINVENT4 repository, which is covered under the Apache License 2.0. You may obtain a copy of the License at:

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Acknowledgements
- RDKit: Open-source cheminformatics
- Contributors and the open-source community

For more details, see the documentation and example notebooks.