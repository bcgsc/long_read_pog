# Dockerfile

FROM jupyter/base-notebook

# Install any additional dependencies you need
# For example, if you use Python notebooks, you might want to install additional packages:
# e.g RUN pip install -r requirements.txt
# e.g RUN pip install pip-install-test

# Change permissions to run jupyter notebook
RUN chmod -R 777 ./

# Set the working directory
WORKDIR /workspace

# Expose the Jupyter notebook port
EXPOSE 8888

# Command to run Jupyter notebook on container startup
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
