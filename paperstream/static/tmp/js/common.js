requirejs.config({
    baseUrl: '/static/tmp/js/lib',
    shim: {
        bootstrap: {
            deps: ['jquery']
        },
        'endless': {
            deps: ['jquery']
        },
        'close_alerts': {
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
        close_alerts: '/static/tmp/js/lib/close-alerts',
        app: '/static/tmp/js/app'
    }
});
