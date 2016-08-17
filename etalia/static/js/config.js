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
        tinymce: {
            exports: 'tinymce'
        }
    },
    map: {
        '*': {
            'backbone/forms': 'lib/backbone/forms-bootstrap'
        },
        'lib/backbone/forms-bootstrap': {
            'backbone-forms': 'lib/backbone/forms'
        }
    },
    packages: [
        {
            name: 'app',
            main: 'app',
            location: 'app'
        },
        {
            name: 'backbone',
            main: 'backbone',
            location: 'lib/backbone'
        }
    ],
    paths: {
        altmetric: 'https://d1bxh8uas1mnw7.cloudfront.net/assets/embed',
        bootstrap: 'lib/bootstrap',
        endless: 'lib/endless-pagination',
        handlebars: 'lib/handlebars',
        hogan: 'lib/hogan',
        jquery: 'lib/jquery',
        'jquery-ui': 'lib/jquery-ui',
        moment: 'lib/moment/moment',
        select2: 'lib/select2/select2',
        tinymce: 'lib/tinymce/tinymce',
        underscore: 'lib/underscore'
    },
    config: {
        moment: {
            noGlobal: true
        }
    }
});
