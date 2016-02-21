var gulp = require('gulp'),
    merge = require('merge-stream'),
    concat = require('gulp-concat'),
    rename = require('gulp-rename'),
    strip = require('gulp-strip-comments'),
    shell = require('gulp-shell'),
    hogan = require('gulp-hogan-compile');



/**
 * (build) and copy libraries
 */
gulp.task('libs', function() {
    /**
     * Require JS
     */
    var requireJs = gulp
        .src('node_modules/requirejs/require.js', {
            base: 'node_modules/requirejs'
        })
        .pipe(strip())
        .pipe(gulp.dest('paperstream/static/js'));

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
        .pipe(gulp.dest('paperstream/static/js/lib'));

    /**
     * jQuery
     */
    var jquery = gulp
        .src('bower_components/jquery/dist/jquery.js', {
            base: 'bower_components/jquery/dist'
        })
        .pipe(strip())
        .pipe(gulp.dest('paperstream/static/js/lib'));

    /**
     * jQuery Ui (custom build)
     * @see http://jqueryui.com/themeroller/?ffDefault=segoe%20ui%2CArial%2Csans-serif&fsDefault=1.1em&fwDefault=bold&cornerRadius=4px&bgColorHeader=%23ffffff&bgTextureHeader=flat&borderColorHeader=%23cbced1&fcHeader=%23313233&iconColorHeader=%23847e71&bgColorContent=%23ffffff&bgTextureContent=flat&borderColorContent=%23cbced1&fcContent=%23313233&iconColorContent=%23808080&bgColorDefault=%2300695c&bgTextureDefault=flat&borderColorDefault=%2300695c&fcDefault=%23ffffff&iconColorDefault=%23eeeeee&bgColorHover=%23004f46&bgTextureHover=flat&borderColorHover=%23004f46&fcHover=%23ffffff&iconColorHover=%23ffffff&bgColorActive=%23fafaf4&bgTextureActive=flat&borderColorActive=%23cbced1&fcActive=%2300695c&iconColorActive=%238DC262&bgColorHighlight=%23cbced1&bgTextureHighlight=glass&borderColorHighlight=%23cbced1&fcHighlight=%23363636&iconColorHighlight=%238DC262&bgColorError=%23ffedad&bgTextureError=highlight_soft&borderColorError=%23e3a345&fcError=%23cd5c0a&iconColorError=%23cd0a0a&bgColorOverlay=%232b2922&bgTextureOverlay=inset_soft&bgImgOpacityOverlay=15&opacityOverlay=90&bgColorShadow=%23cccccc&bgTextureShadow=highlight_hard&bgImgOpacityShadow=95&opacityShadow=20&thicknessShadow=12px&offsetTopShadow=-12px&offsetLeftShadow=-12px&cornerRadiusShadow=10px&bgImgOpacityHeader=0&bgImgOpacityContent=0&bgImgOpacityDefault=0&bgImgOpacityHover=0&bgImgOpacityActive=0&bgImgOpacityHighlight=55&bgImgOpacityError=95
     * core.js, widget.js, mouse.js, position.js, slider.js
     */
    var jqueryUi = gulp
        .src([
            'bower_components/jquery-ui/ui/core.js',
            'bower_components/jquery-ui/ui/widget.js',
            'bower_components/jquery-ui/ui/mouse.js',
            'bower_components/jquery-ui/ui/position.js',
            'bower_components/jquery-ui/ui/slider.js'
        ])
        .pipe(concat('jquery-ui.js'))
        .pipe(strip())
        .pipe(gulp.dest('paperstream/static/js/lib'));

    /**
     * jQuery mouse wheel
     */
    var jqueryMouseWheel = gulp
        .src('bower_components/jquery-mousewheel/jquery.mousewheel.js', {
            base: 'bower_components/jquery-mousewheel'
        })
        .pipe(strip())
        .pipe(gulp.dest('paperstream/static/js/lib'));

    /**
     * Bootstrap (custom build)
     */
    var bootstrap = gulp
        .src([
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/affix.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/alert.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/collapse.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/modal.js',
            'bower_components/bootstrap-sass/assets/javascripts/bootstrap/transition.js'
        ])
        .pipe(concat('bootstrap.js'))
        .pipe(strip())
        .pipe(gulp.dest('paperstream/static/js/lib'));

    return merge(requireJs, ie9, jquery, jqueryUi, jqueryMouseWheel, bootstrap);
});



/**
 * Bootstrap SASS
 */
// TODO Bootstrap sass
// http://david-barreto.com/working-with-sass-bootstrap-and-gulp/



/**
 * Mustache templates compilation
 */
gulp.task('hogan', function() {
    return gulp
        .src('paperstream/static/templates/*.mustache')
        .pipe(hogan('templates.js', {
            wrapper: 'amd',
            hoganModule: 'hogan'
        }))
        .pipe(gulp.dest('paperstream/static/js/app'));
});



/**
 * RequireJs optimisation
 */
gulp.task('rjs', shell.task([
    'node_modules/requirejs/bin/r.js -o paperstream/static/js/build.js'
]));


/**
 * Tasks
 */
gulp.task('default', ['hogan', 'rjs']);
