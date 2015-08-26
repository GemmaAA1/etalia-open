#!/usr/bin/env bash

# Check if virtual env is active
if [ -z "$VIRTUAL_ENV" ]; then
    echo 'VIRTUAL_ENV not set'
    exit 1
fi

# export environment variable to virtulalenv postactivate
file=$VIRTUAL_ENV/bin/postactivate
rm $file
echo "#!/usr/bin/env bash" >> $file
echo "export DJANGO_LOG_LEVEL='DEBUG'" >> $file
echo "export DJANGO_DEBUG=True" >> $file
#echo "export DJANGO_EMAIL_BACKEND=''" >> $file
echo "export PAP_SECRET_KEY='_kjjr3)+tdlfbcw7uu&oue+*50+hbv9gsd-yx35^*%n\$5ugp-s'" >> $file
echo "export ELSEVIER_API_KEY='293e288c325d7765b7c22f5195175351'" >> $file
echo "export PUBMED_EMAIL='nicolas.pannetier@gmail.com'" >> $file
source $file

## cleaning up symmetrically
file=$VIRTUAL_ENV/bin/predeactivate
echo "unset DJANGO_LOG_LEVEL" >> $file
echo "unset DJANGO_DEBUG" >> $file
#echo "unset DJANGO_EMAIL_BACKEND" >> $file
echo "unset PAP_SECRET_KEY" >> $file
echo "unset ELSEVIER_API_KEY" >> $file
echo "unset PUBMED_EMAIL" >> $file