FROM python:3.4

ENV PYTHONUNBUFFERED 1

COPY requirements /requirements

RUN pip install -r /requirements/development.txt
RUN adduser --disabled-password --gecos '' myuser
RUN mkdir /log
RUN mkdir /data
RUN mkdir /code
RUN chown -R myuser:myuser /log
RUN chown -R myuser:myuser /data
RUN chown -R myuser:myuser /root

WORKDIR /code
