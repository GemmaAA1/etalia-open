var gulp = require('gulp'),
    hogan = require('gulp-hogan-compile');

// TODO copy libs

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
