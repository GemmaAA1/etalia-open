FROM python:3.4
ENV PYTHONUNBUFFERED 1
RUN mkdir /code
WORKDIR /code
COPY requirements /code/requirements
RUN pip install -r requirements/development.txt
ADD . /code/
RUN adduser --disabled-password --gecos '' myuser
RUN mkdir /log
RUN mkdir /data
RUN chown -R myuser:myuser /log
RUN chown -R myuser:myuser /data
