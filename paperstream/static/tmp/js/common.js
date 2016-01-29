requirejs.config({
    baseUrl: '/static/tmp/js/lib',
    shim: {
        bootstrap: {
            deps: ['jquery']
        },
        'endless': {
            deps: ['jquery']
        },
        'jqtip': {
            deps: ['jquery']
        },
        'jquery-ui': {
            deps: ['jquery']
        }
    },
    paths: {
        endless: '/static/tmp/js/lib/endless-pagination',
        jqtip: '/static/tmp/js/lib/jquery.qtip',
        app: '/static/tmp/js/app'
    }
});
