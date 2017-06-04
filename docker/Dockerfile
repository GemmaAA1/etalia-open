FROM python:3.4

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

CMD python manage.py runserver 0.0.0.0:8000