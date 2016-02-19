requirejs.config({
    baseUrl: '/static/js/lib',
    shim: {
        bootstrap: {
            deps: ['jquery']
        },
        'endless': {
            deps: ['jquery']
        },
        'close-alerts': {
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
        endless: '/static/js/lib/endless-pagination',
        app: '/static/js/app'
    }
});
