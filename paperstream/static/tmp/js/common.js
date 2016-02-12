requirejs.config({
    baseUrl: '/static/tmp/js/lib',
    shim: {
        bootstrap: {
            deps: ['jquery']
        },
        'endless': {
            deps: ['jquery']
        },
        'jquery-ui': {
            deps: ['jquery']
        },
        'jquery.mousewheel': {
            deps: ['jquery']
        }
    },
    paths: {
        endless: '/static/tmp/js/lib/endless-pagination',
        app: '/static/tmp/js/app'
    }
});
