#!/bin/bash
spark-submit --name test --master local --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=hdfs://node001:8020/user/digest/logs GVR.py
sshpass -p 'crayon99' scp *.txt yleborgn@bridgeiris.ulb.ac.be:/home/yleborgn/digest/users/guest/analyses
