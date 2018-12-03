FROM mjstealey/irods-icommands:4.1.10

RUN mkdir /cvmfs /work
WORKDIR /work

# Install some prerequisites.
RUN apt-get update \
    && apt-get install -y lsb wget apt-transport-https python2.7 python-requests

# Install the wrapper script and the script to upload the output files.
ADD wrapper /usr/bin/wrapper
ADD upload-files /usr/bin/upload-files

# Make the wrapper script the default command.
CMD ["wrapper"]
