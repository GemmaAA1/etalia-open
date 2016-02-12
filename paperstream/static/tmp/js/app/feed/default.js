define(['jquery', 'app/api', 'app/ui/list', 'bootstrap'], function($, Api) {

    $(function() {

        $('body')
            .on('etalia.publication.pin', function(e, data) {
                if (data.hasOwnProperty('is_pinned')) {
                    $('.thumb[data-id=' + data.id + ']')
                        .find('.thumb-pin')
                        .toggleClass('active', data['is_pinned']);
                }
            })
            .on('etalia.publication.ban', function(e, data) {
                if (data.hasOwnProperty('is_banned') && data['is_banned']) {
                    /*$('.thumb[data-id=' + data.id + ']')
                     .find('.thumb-ban')
                     .toggleClass('active', data['success']);*/

                    $('.thumb[data-id=' + data.id + ']').remove();
                }
            })
            .on('etalia.publication.trash', function(e, data) {
                if (data.hasOwnProperty('success') && data['success']) {
                    $('.thumb[data-id=' + data.id + ']').remove();
                }
            });
    });
});
