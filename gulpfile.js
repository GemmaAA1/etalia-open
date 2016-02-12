var gulp = require('gulp'),
    hogan = require('gulp-hogan-compile');

// Copy libs
gulp.task('libs', function() {
    gulp.src('node_modules/jquery-mousewheel/jquery.mousewheel.js', {base: 'node_modules/jquery-mousewheel'})
        .pipe(gulp.dest('paperstream/static/tmp/js/lib'));
});

// TODO Bootstrap sass

// Mustache templates compilation
gulp.task('hogan', function() {
    gulp.src('paperstream/static/tmp/templates/*.mustache')
        .pipe(hogan('templates.js', {
            wrapper: 'amd',
            hoganModule: 'hogan'
        }))
        .pipe(gulp.dest('paperstream/static/tmp/js/app'));
});

gulp.task('default', ['hogan']);
