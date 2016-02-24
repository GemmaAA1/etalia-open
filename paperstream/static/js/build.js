//noinspection BadExpressionStatementJS
({
    // Relative to this file's location
    baseUrl: '.',
    // If you got a main config file, this is the palce for it
    mainConfigFile: 'config.js',
    // The destination directory
    dir: '../compiled/js',
    // If set to true, any files that were combined into a
    // build bundle will be removed from the output folder.
    removeCombined: true,
    // Finds require() dependencies inside a require() or define call. By default
    // this value is false, because those resources should be considered dynamic/runtime
    // calls. However, for some optimization scenarios, it is desirable to
    // include them in the build.
    //findNestedDependencies: false,
    // For test purposes, the output is not minified.
    // In the real world this would be set to uglify, uglify2, or to closure
    optimize: 'uglify2',

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
        }
    ]
})
