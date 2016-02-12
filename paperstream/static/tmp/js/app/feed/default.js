define(['jquery', 'app/api', 'app/ui/list', 'bootstrap'], function($) {

    $(function() {

        $('body')
            .on('etalia.publication.pin', function(e, result) {
                $('.thumb[data-id=' + result.getId() + ']')
                    .find('.thumb-pin')
                    .toggleClass('active', result.isPinned());
            })
            .on('etalia.publication.ban', function(e, result) {
                if (result.isBanned()) {
                    $('.thumb[data-id=' + result.getId() + ']').remove();
                }
            })
            .on('etalia.publication.trash', function(e, result) {
                if (result.isTrashed()) {
                    $('.thumb[data-id=' + result.getId() + ']').remove();
                }
            });
    });
});
