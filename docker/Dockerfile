FROM python:3.4.3

RUN apt-get -y update
RUN apt-get -y install libc-dev
RUN apt-get -y install build-essential
RUN pip install -U pip

COPY requirements /requirements

RUN pip install -r /requirements/development.txt \
    && adduser --disabled-password --gecos '' myuser \
    && mkdir /log \
    && mkdir /data \
    && mkdir /code \
    && chown -R myuser:myuser /log \
    && chown -R myuser:myuser /data \
    && chown -R myuser:myuser /root

WORKDIR /code

RUN apt-get autoremove

CMD python manage.py runserver 0.0.0.0:8000