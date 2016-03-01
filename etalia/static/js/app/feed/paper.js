define([
    'jquery',
    'app/ui/paper',
    'app/ui/layout',
    'bootstrap'
], function($, Paper) {

    $(function() {

        var paper = new Paper();
        paper.init().update();

        $('body')
            .on('etalia.publication.ban', function(e, result) {
                if (result.isBanned()) {
                    window.location.href = '/feed/stream';
                }
            });
    });
});
