FROM python:3.4.3

RUN apt-get -y update \
 && apt-get -y install libc-dev build-essential \
 && pip install -U pip

COPY requirements /requirements

RUN pip install -r /requirements/development.txt \
    && adduser --disabled-password --gecos '' myuser \
    && mkdir /log \
    && mkdir /data \
    && mkdir /code \
    && chown -R myuser:myuser /log \
    && chown -R myuser:myuser /data \
    && chown -R myuser:myuser /root

RUN apt-get autoremove -y --purge \
    && apt-get clean -y

WORKDIR /code

CMD python manage.py runserver 0.0.0.0:8000
