FROM mjstealey/irods-icommands:4.1.10

RUN mkdir /cvmfs /work
WORKDIR /work

# Install some prerequisites.
RUN apt-get update \
    && apt-get install -y jq lsb wget curl apt-transport-https python2.7 python-requests

# Install the wrapper script and the script to upload the output files.
COPY wrapper /usr/bin/wrapper
COPY send-job-status-update /usr/bin/send-job-status-update
COPY ticketParser.sh /usr/bin/ticketParser.sh
COPY evalTicket.sh /usr/bin/evalTicket.sh

# Make the wrapper script the default command.
CMD ["wrapper"]
