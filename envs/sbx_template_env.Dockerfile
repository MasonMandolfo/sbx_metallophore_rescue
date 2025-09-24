FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_metallophore_rescue_env

COPY envs/sbx_metallophore_rescue_env.yml ./

# Install environment
RUN conda env create --file sbx_metallophore_rescue_env.yml --name sbx_metallophore_rescue

ENV PATH="/opt/conda/envs/sbx_metallophore_rescue/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_metallophore_rescue", "/bin/bash", "-c"]

# Run
CMD "bash"