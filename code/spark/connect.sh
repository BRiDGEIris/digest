#!/bin/bash
scp -P 23 /srv/shiny-server/digest/users/guest/jobsArguments.conf digest@hadoop.mlg.ulb.ac.be:/home/digest
ssh -p 23 -i /home/digest/.ssh/id_rsa digest@hadoop.mlg.ulb.ac.be 'bash -s' < /srv/shiny-server/digest/run.sh

