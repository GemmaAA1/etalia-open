define(['jquery', 'app/ui/paper', 'app/ui/layout', 'bootstrap'], function($, Paper) {

    $(function() {

        var paper = new Paper();
        paper.init();

        $('body')
            .on('etalia.publication.ban', function(e, result) {
                if (result.isBanned()) {
                    // TODO redirect ?
                }
            });
    });
});
