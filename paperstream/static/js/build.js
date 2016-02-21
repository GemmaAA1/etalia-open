//noinspection BadExpressionStatementJS
({
    // Relative to this file's location
    baseUrl: '.',
    // If you got a main config file, this is the palce for it
    mainConfigFile: 'common.js',
    // The destination directory
    dir: '../../static-dist/js',
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
            name: 'app/util/templates'
        },
        {
            name: 'app/util/utils',
            exclude: ['jquery']
        },
        {
            name: 'app/ui/controls',
            exclude: ['jquery', 'app/util/templates', 'app/util/utils']
        },
        {
            name: 'app/ui/layout',
            exclude: ['jquery', 'app/util/utils']
        },
        {
            name: 'app/ui/paper',
            exclude: ['jquery', 'app/api', 'app/util/utils', 'app/ui/controls', 'app/ui/list']
        }
    ]
})
