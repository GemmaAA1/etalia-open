FROM python:3.4


ENV PYTHONUNBUFFERED 1


RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        postgresql-client \
    && rm -rf /var/lib/apt/lists/*


EXPOSE 8000


RUN mkdir /etalia
WORKDIR /etalia


# Copy the {host_project_root}/requirements folder into the container /requirements folder
COPY requirements /requirements

RUN pip3 install -r /requirements/development.txt
RUN pip install --upgrade numpy


COPY ./docker/engine/entrypoint.sh /entrypoint.sh

RUN chmod +x /entrypoint.sh
