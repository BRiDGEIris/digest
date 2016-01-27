FROM rocker/shiny

#Install Java and a few tools
RUN apt-get update 
RUN apt-get -y install libmysqlclient-dev libxml2-dev libssl-dev git nano
RUN apt-get install -y openjdk-7-jdk

#Install R packages required for DiGeST
#Basic
RUN R -e "install.packages(c('jsonlite','plyr', 'ggplot2','devtools'), repos='http://cran.rstudio.com/')"
#HTML widgets
RUN R -e "install.packages(c('htmlwidgets','DT','shinyBS', 'shinyjs'), repos='http://cran.rstudio.com/')"

ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre
ENV PATH $PATH:$JAVA_HOME/bin
RUN R CMD javareconf
#Database connectors
RUN R -e "install.packages(c('RMySQL','RSQLite','RJDBC','RCurl'), repos='http://cran.rstudio.com/')"

# queryBuildR and rpivotTable need to be installed from GitHub
WORKDIR /home/docker
RUN git clone https://github.com/Yannael/queryBuildR.git
RUN R CMD INSTALL queryBuildR
RUN git clone https://github.com/smartinsightsfromdata/rpivotTable.git
RUN R CMD INSTALL rpivotTable

RUN apt-get update
RUN apt-get install -y openvpn
ADD yleborgn-vpn-cc.vub.ac.be.tar.gz /home/docker/openvpn
RUN mv /home/docker/openvpn/* /etc/openvpn

ADD connectHighlander.R /srv/shiny-server/connectHighlander.R

RUN apt-get install -y net-tools whois
#RUN mkdir /root/.ssh
#RUN ssh-keygen -t rsa -N '' -f /root/.ssh/id_rsa

#RUN useradd -m digest -p `mkpasswd digest_1239`
#USER digest
#WORKDIR /home/digest
#RUN mkdir .ssh
RUN mkdir /home/shiny/.ssh
ADD id_rsa /home/shiny/.ssh/id_rsa
ADD known_hosts /home/shiny/.ssh/known_hosts

RUN echo shiny:shiny1239 | /usr/sbin/chpasswd

RUN apt-get install -y openssh-server 

ADD connect.sh /srv/shiny-server/digest/connect.sh
ADD run.sh /srv/shiny-server/digest/run.sh

ADD startupscript.sh /usr/bin/startupscript.sh

RUN R -e "devtools::install_github('rstudio/DT')"


#Expose ports for Spark UI and history server
EXPOSE 1194
EXPOSE 21050
EXPOSE 22

CMD ["shiny-server"]

