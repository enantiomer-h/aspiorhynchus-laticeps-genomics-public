# NOTE FOR EXTERNAL USERS
# ------------------------------------------------------------------------------
# `huang/comparative-genomics:20250828` is an INTERNAL image and is not published
# to a public registry, so this Dockerfile cannot be built as-is outside the
# authors' environment. It is retained only to document the original container.
#
# The supported, reproducible setup is the Conda environments in `env/*.yml`
# (see the repository README, "Environment setup"). Create them with, e.g.:
#     conda env create -f env/env_biopython.yml
# and select the matching environment per notebook. Any base Jupyter image with
# these Conda environments installed reproduces the analysis.
# ------------------------------------------------------------------------------
FROM huang/comparative-genomics:20250828

USER root

CMD echo 'Example for Starting Jupyter Notebook on port 30000: jupyter notebook --ip="*" --port=30000 --no-browser --allow-root' && \
    jupyter notebook \
    --ip='*' --port=30000 \
    --no-browser --allow-root
