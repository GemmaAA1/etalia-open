requirejs.config({
    baseUrl: '/static/tmp/js/lib',
    shim: {
        bootstrap: {
            deps: ['jquery']
        }
    },
    paths: {
        app: '/static/tmp/js/app'
    }
});
