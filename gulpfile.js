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
 * (build) and copy libraries
 */
gulp.task('libraries', function() {
    /**
     * Require JS
     */
    var requireJs = gulp
        .src('node_modules/requirejs/require.js', {
            base: 'node_modules/requirejs'
        })
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js'));

    /**
     * ie9 (html5shiv + respond)
     */
    var ie9 = gulp
        .src([
            'bower_components/html5shiv/dist/html5shiv.js',
            'bower_components/respond/dest/respond.src.js'
        ])
        .pipe(concat('ie9.js'))
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js/lib'));

    /**
     * jQuery
     */
    var jquery = gulp
        .src('bower_components/jquery/dist/jquery.js', {
            base: 'bower_components/jquery/dist'
        })
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js/lib'));

    /**
     * jQuery Ui (custom build)
     * @see http://jqueryui.com/themeroller/?ffDefault=segoe%20ui%2CArial%2Csans-serif&fsDefault=1.1em&fwDefault=bold&cornerRadius=4px&bgColorHeader=%23ffffff&bgTextureHeader=flat&borderColorHeader=%23cbced1&fcHeader=%23313233&iconColorHeader=%23847e71&bgColorContent=%23ffffff&bgTextureContent=flat&borderColorContent=%23cbced1&fcContent=%23313233&iconColorContent=%23808080&bgColorDefault=%2300695c&bgTextureDefault=flat&borderColorDefault=%2300695c&fcDefault=%23ffffff&iconColorDefault=%23eeeeee&bgColorHover=%23004f46&bgTextureHover=flat&borderColorHover=%23004f46&fcHover=%23ffffff&iconColorHover=%23ffffff&bgColorActive=%23fafaf4&bgTextureActive=flat&borderColorActive=%23cbced1&fcActive=%2300695c&iconColorActive=%238DC262&bgColorHighlight=%23cbced1&bgTextureHighlight=glass&borderColorHighlight=%23cbced1&fcHighlight=%23363636&iconColorHighlight=%238DC262&bgColorError=%23ffedad&bgTextureError=highlight_soft&borderColorError=%23e3a345&fcError=%23cd5c0a&iconColorError=%23cd0a0a&bgColorOverlay=%232b2922&bgTextureOverlay=inset_soft&bgImgOpacityOverlay=15&opacityOverlay=90&bgColorShadow=%23cccccc&bgTextureShadow=highlight_hard&bgImgOpacityShadow=95&opacityShadow=20&thicknessShadow=12px&offsetTopShadow=-12px&offsetLeftShadow=-12px&cornerRadiusShadow=10px&bgImgOpacityHeader=0&bgImgOpacityContent=0&bgImgOpacityDefault=0&bgImgOpacityHover=0&bgImgOpacityActive=0&bgImgOpacityHighlight=55&bgImgOpacityError=95
     * core.js, widget.js, mouse.js, slider.js
     */
    var jqueryUi = gulp
        .src([
            'bower_components/jquery-ui/ui/core.js',
            'bower_components/jquery-ui/ui/widget.js',
            'bower_components/jquery-ui/ui/mouse.js',
            'bower_components/jquery-ui/ui/slider.js'
        ], {
            base: 'bower_components/jquery-ui/ui'
        })
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js/lib/jquery-ui'));

    /**
     * Bootstrap (custom build)
     */
    var bootstrap = gulp
        .src([
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/affix.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/alert.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/collapse.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/modal.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/tooltip.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/transition.js'
        ])
        .pipe(concat('bootstrap.js'))
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js/lib'));

    /**
     * Hogan.js
     */
    var hogan = gulp
        .src('node_modules/hogan.js/dist/hogan-3.0.2.amd.js', {
            base: 'node_modules/hogan.js/dist'
        })
        .pipe(rename('hogan.js'))
        .pipe(strip())
        .pipe(gulp.dest(config.src + '/js/lib'));

    return merge(requireJs, ie9, jquery, jqueryUi, bootstrap, hogan);
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
            config.src + '/css/app/root.css'
        ])
        .pipe(concat('main.css'))
        .pipe(replace(/..\/..\/([a-z]+)/g, '../$1')) // Fix bootstrap and jquery-ui fonts/img paths
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
            config.src + '/css/app/thread.css',
            config.src + '/css/app/thread-content.css'
        ])
        .pipe(concat('elements.css'))
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    var user = gulp
        .src(config.src + '/css/app/user.css')
        .pipe(minify({compatibility: 'ie8'}))
        .pipe(gulp.dest(config.dest + '/css'));

    return merge(main, page, landing, elements, user);
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
