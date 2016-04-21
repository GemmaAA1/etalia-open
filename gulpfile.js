var gulp = require('gulp'),
    clean = require('gulp-clean'),
    merge = require('merge-stream'),
    concat = require('gulp-concat'),
    rename = require('gulp-rename'),
    replace = require('gulp-replace'),
    strip = require('gulp-strip-comments'),
    shell = require('gulp-shell'),
    hogan = require('gulp-hogan-compile'),
    minify = require('gulp-minify-css'),
    runSequence = require('run-sequence');

var config = {
    src: 'etalia/static',
    dest: 'etalia/static/compiled'
};



/**
 * Clean compiled directory.
 */
gulp.task('clean-libraries', function () {
    return gulp
        .src(config.src + '/js/lib', {read: false})
        .pipe(clean());
});

/**
 * (build) and copy libraries
 */
gulp.task('libraries', ['clean-libraries'], function() {
    var libs = [
        {
            src: 'node_modules/requirejs/require.js',
            base: 'node_modules/requirejs',
            dest: '/js'
        },
        {
            src: 'bower_components/text/text.js',
            base: 'bower_components/text',
            dest: '/js'
        },
        {
            src: [
                'bower_components/html5shiv/dist/html5shiv.js',
                'bower_components/respond/dest/respond.src.js'
            ],
            concat: 'ie9.js',
            base: '',
            dest: '/js/lib'
        },
        {
            src: 'bower_components/underscore/underscore.js',
            base: 'bower_components/underscore',
            dest: '/js/lib'
        },
        {
            src: 'bower_components/backbone/backbone.js',
            base: 'bower_components/backbone',
            dest: '/js/lib/backbone'
        },
        {
            src: 'bower_components/backbone-relational/backbone-relational.js',
            base: 'bower_components/backbone-relational',
            rename: 'relational.js',
            dest: '/js/lib/backbone'
        },
        {
            src: 'bower_components/backbone.paginator/lib/backbone.paginator.js',
            base: 'bower_components/backbone.paginator',
            rename: 'paginator.js',
            dest: '/js/lib/backbone'
        },
        {
            src: 'bower_components/backbone-forms/distribution.amd/backbone-forms.js',
            base: 'bower_components/backbone-forms/distribution.amd',
            rename: 'forms.js',
            dest: '/js/lib/backbone'
        },
        {
            src: 'bower_components/backbone-forms/distribution.amd/templates/bootstrap3.js',
            base: 'bower_components/backbone-forms/distribution.amd/templates',
            rename: 'forms-bootstrap.js',
            dest: '/js/lib/backbone'
        },
        {
            src: 'bower_components/backbone-forms/distribution.amd/templates/bootstrap3.css',
            base: 'bower_components/backbone-forms/distribution.amd/templates',
            rename: 'backbone-forms-bootstrap.css',
            dest: '/css/lib'
        },
        {
            src: 'bower_components/handlebars/handlebars.js',
            base: 'bower_components/handlebars',
            dest: '/js/lib'
        },
        {
            src: 'bower_components/jquery/dist/jquery.js',
            base: 'bower_components/jquery/dist',
            dest: '/js/lib'
        },
        {
            // @see http://jqueryui.com/themeroller/?ffDefault=segoe%20ui%2CArial%2Csans-serif&fsDefault=1.1em&fwDefault=bold&cornerRadius=4px&bgColorHeader=%23ffffff&bgTextureHeader=flat&borderColorHeader=%23cbced1&fcHeader=%23313233&iconColorHeader=%23847e71&bgColorContent=%23ffffff&bgTextureContent=flat&borderColorContent=%23cbced1&fcContent=%23313233&iconColorContent=%23808080&bgColorDefault=%2300695c&bgTextureDefault=flat&borderColorDefault=%2300695c&fcDefault=%23ffffff&iconColorDefault=%23eeeeee&bgColorHover=%23004f46&bgTextureHover=flat&borderColorHover=%23004f46&fcHover=%23ffffff&iconColorHover=%23ffffff&bgColorActive=%23fafaf4&bgTextureActive=flat&borderColorActive=%23cbced1&fcActive=%2300695c&iconColorActive=%238DC262&bgColorHighlight=%23cbced1&bgTextureHighlight=glass&borderColorHighlight=%23cbced1&fcHighlight=%23363636&iconColorHighlight=%238DC262&bgColorError=%23ffedad&bgTextureError=highlight_soft&borderColorError=%23e3a345&fcError=%23cd5c0a&iconColorError=%23cd0a0a&bgColorOverlay=%232b2922&bgTextureOverlay=inset_soft&bgImgOpacityOverlay=15&opacityOverlay=90&bgColorShadow=%23cccccc&bgTextureShadow=highlight_hard&bgImgOpacityShadow=95&opacityShadow=20&thicknessShadow=12px&offsetTopShadow=-12px&offsetLeftShadow=-12px&cornerRadiusShadow=10px&bgImgOpacityHeader=0&bgImgOpacityContent=0&bgImgOpacityDefault=0&bgImgOpacityHover=0&bgImgOpacityActive=0&bgImgOpacityHighlight=55&bgImgOpacityError=95
            src: [
                'bower_components/jquery-ui/ui/core.js',
                'bower_components/jquery-ui/ui/widget.js',
                'bower_components/jquery-ui/ui/mouse.js',
                'bower_components/jquery-ui/ui/slider.js'
            ],
            base: 'bower_components/jquery-ui/ui',
            dest: '/js/lib/jquery-ui'
        },
        {
            src: [
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/affix.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/alert.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/button.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/collapse.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/modal.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/tooltip.js',
                'bower_components/bootstrap-sass/assets/javascripts/bootstrap/transition.js'
            ],
            base: 'bower_components/bootstrap-sass/assets/javascripts/bootstrap',
            concat: 'bootstrap.js',
            dest: '/js/lib'
        },
        {
            src: [
                'bower_components/tinymce/plugins/lists/plugin.js',
                'bower_components/tinymce/plugins/paste/plugin.js',
                'bower_components/tinymce/plugins/link/plugin.js',
                'bower_components/tinymce/plugins/searchreplace/plugin.js',
                'bower_components/tinymce/plugins/hr/plugin.js',
                //'bower_components/tinymce/skins/**/*',
                'bower_components/tinymce/themes/**/*',
                'bower_components/tinymce/tinymce.js'
            ],
            base: 'bower_components/tinymce',
            dest: '/js/lib/tinymce'
        },
        {
            src: [
                'bower_components/moment/moment.js'
                //'bower_components/moment/locale/*'
            ],
            base: 'bower_components/moment',
            dest: '/js/lib/moment'
        },
        {
            src: 'node_modules/hogan.js/dist/hogan-3.0.2.amd.js',
            base: 'node_modules/hogan.js/dist',
            rename: 'hogan.js',
            dest: '/js/lib'
        }
    ];

    var merged = merge();
    for (var i in libs) {
        var lib = libs[i],
            stream;
        if (lib.base) {
            stream = gulp.src(lib.src, {base: lib.base});
        } else {
            stream = gulp.src(lib.src);
        }
        if (lib.concat) {
            stream = stream.pipe(concat(lib.concat));
        }
        if (lib.rename) {
            stream = stream.pipe(rename(lib.rename));
        }
        stream = stream
            //.pipe(strip()) // @TODO bug with css
            .pipe(gulp.dest(config.src + lib.dest));
        merged.add(stream);
    }

    return merged;
});



/**
 * Mustache templates compilation
 */
gulp.task('templates', function() {
    return gulp
        .src(config.src + '/templates/*.mustache')
        .pipe(hogan('templates.js', {
            wrapper: 'amd',
            hoganModule: 'hogan'
        }))
        .pipe(gulp.dest(config.src + '/js/app/util'));
});



/**
 * RequireJs optimisation
 */
gulp.task('scripts', shell.task([
    'node_modules/requirejs/bin/r.js -o etalia/static/js/build.js'
]));



/**
 * Bootstrap SASS
 */
// TODO Bootstrap sass
// http://david-barreto.com/working-with-sass-bootstrap-and-gulp/



/**
 * Minify css
 */
gulp.task('styles', function() {
    var main = gulp
        .src([
            config.src + '/css/*.css',
            config.src + '/css/lib/*.css',
            config.src + '/css/app/content.css',
            config.src + '/css/app/root.css'
        ])
        .pipe(concat('main.css'))
        .pipe(replace(/..\/..\/([a-z]+)/g, '../$1')) // Fix bootstrap and jquery-ui fonts/img paths
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var content = gulp
        .src(config.src + '/css/app/content.css')
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var page = gulp
        .src(config.src + '/css/app/page.css')
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var landing = gulp
        .src(config.src + '/css/app/landing.css')
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var elements = gulp
        .src([
            config.src + '/css/app/list.css',
            config.src + '/css/app/card.css',
            config.src + '/css/app/detail.css',
            config.src + '/css/app/thread.css'
        ])
        .pipe(concat('elements.css'))
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var user = gulp
        .src(config.src + '/css/app/user.css')
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    return merge(main, content, page, landing, elements, user);
});



/**
 * Images and fonts
 */
gulp.task('images', function() {
    return gulp
        .src(config.src + '/img/**')
        .pipe(gulp.dest(config.dest + '/img'));
});
gulp.task('fonts', function() {
    return gulp
        .src(config.src + '/fonts/**')
        .pipe(gulp.dest(config.dest + '/fonts'));
});



/**
 * Clean compiled directory.
 */
gulp.task('clean', function () {
    return gulp
        .src(config.dest, {read: false})
        .pipe(clean());
});



/**
 * Tasks
 */
gulp.task('build', function() {
    runSequence(
        'clean',
        ['libraries', 'templates'],
        ['scripts', 'styles', 'images', 'fonts']
    );
});
gulp.task('default', ['build']);
