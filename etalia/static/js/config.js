requirejs.config({
    baseUrl: '/static/js',
    shim: {
        bootstrap: {
            deps: ['jquery']
        },
        'endless': {
            deps: ['jquery']
        },
        'jquery-ui': {
            deps: ['jquery']
        }
    },
    paths: {
        bootstrap: './lib/bootstrap',
        endless: './lib/endless-pagination',
        hogan: './lib/hogan',
        jquery: './lib/jquery',
        'jquery-ui': './lib/jquery-ui',
        altmetric: 'https://d1bxh8uas1mnw7.cloudfront.net/assets/embed'
    }
});