#!/bin/bash
#/etc/init.d/openvpn restart
bash << EOF
shiny-server&
ssh -nNT -p 23 -L 21050:127.0.0.1:21050 digest@hadoop.mlg.ulb.ac.be&
EOF

