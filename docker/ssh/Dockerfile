FROM etalia/dev

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install openssh-server sudo
ADD docker/ssh/set_root_pw.sh /set_root_pw.sh
ADD docker/ssh/run.sh /run.sh
ADD docker/ssh/ssh_keys/myuser-sshd.pub /myuser-sshd.pub
RUN chmod +x /*.sh
RUN mkdir -p /var/run/sshd && sed -i "s/UsePrivilegeSeparation.*/UsePrivilegeSeparation no/g" /etc/ssh/sshd_config \
  && touch /root/.Xauthority \
  && true

## Set a default user. Available via runtime flag `--user docker`
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory, but also be able to sudo
#RUN useradd myuser \
#        && passwd -d myuser \
#        && mkdir /home/myuser \
#        && chown myuser:myuser /home/myuser \
#        && addgroup myuser staff \
#        && addgroup myuser sudo \
#        && true

EXPOSE 22
CMD ["/run.sh"]