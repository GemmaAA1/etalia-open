//noinspection BadExpressionStatementJS
({
    baseUrl: '.',
    mainConfigFile: 'config.js',
    dir: '../compiled/js',
    removeCombined: true,
    optimize: 'uglify2',

    paths: {
        altmetric: 'empty:',
        userlike: 'empty:'
    },

    modules: [
        {
            name: 'app/default',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/pages',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/feed/default',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/feed/paper',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/user/profile',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/user/settings',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/user/signup',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app/user/library',
            exclude: ['jquery', 'bootstrap']
        },
        {
            name: 'app',
            exclude: ['require', 'jquery', 'underscore', 'backbone', 'handlebars', 'moment', 'bootstrap', 'select2']
        },
        {
            name: 'app/feed/list',
            exclude: ['app/app', 'select2']
        },
        {
            name: 'app/paper/list',
            exclude: ['app/app', 'select2']
        },
        {
            name: 'app/thread/list',
            exclude: ['app/app', 'select2']
        }
    ]
})
