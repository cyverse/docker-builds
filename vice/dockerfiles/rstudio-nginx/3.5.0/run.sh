#!/bin/sh

gomplate -f /nginx.conf.tmpl -o /etc/nginx/nginx.conf
supervisord -c /etc/supervisor/supervisord.conf -n
