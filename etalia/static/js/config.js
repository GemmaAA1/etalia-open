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
        },
        'backbone-modal': {
            deps: ['backbone']
        }
    },
    map: {
        '*': {
            'backbone-forms': 'lib/backbone-forms-bootstrap'
        },
        'lib/backbone-forms-bootstrap': {
            'backbone-forms': 'lib/extend/backbone-forms'
        }
    },
    paths: {
        altmetric: 'https://d1bxh8uas1mnw7.cloudfront.net/assets/embed',
        backbone: 'lib/extend/backbone',
        'backbone-modal': 'lib/backbone-modal',
        'backbone-relational': 'lib/backbone-relational',
        bootstrap: 'lib/bootstrap',
        endless: 'lib/endless-pagination',
        handlebars: 'lib/handlebars',
        hogan: 'lib/hogan',
        jquery: 'lib/jquery',
        'jquery-ui': 'lib/jquery-ui',
        underscore: 'lib/underscore'
    }
});
