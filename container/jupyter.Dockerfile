# Dockerfile

FROM jupyter/base-notebook



# must reset to user root to install packages
USER root

# Install Arial font and co:
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections
RUN apt-get update && apt-get install -y font-manager \
 msttcorefonts && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*

# set the user back to original setting
USER $NB_UID

# Install any additional dependencies you need
# For example, if you use Python notebooks, you might want to install additional packages:
# e.g RUN pip install -r requirements.txt
# e.g RUN pip install pip-install-test
RUN pip install pandas 
RUN pip install scipy
RUN pip install datashader
RUN pip install matplotlib
RUN pip install seaborn
RUN pip install numpy
RUN pip install statsmodels

# Change permissions to run jupyter notebook
RUN chmod -R 777 ./

# Set the working directory
WORKDIR /workspace

# Expose the Jupyter notebook port
EXPOSE 8888

# Command to run Jupyter notebook on container startup
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
