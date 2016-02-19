requirejs.config({
    baseUrl: '/static/js/lib',
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
        endless: '/static/js/lib/endless-pagination',
        close_alerts: '/static/js/lib/close-alerts',
        app: '/static/js/app'
    }
});
