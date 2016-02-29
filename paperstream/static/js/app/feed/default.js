define([
    'jquery',
    'app/ui/layout',
    'app/ui/controls',
    'app/ui/list',
    'app/ui/detail',
    'endless',
    'bootstrap'
], function($, layout, controls, List, Detail) {

    var $body,
        list = new List(),
        detail = new Detail();

    $(function() {

        $body = $('body');

        list.init();
        detail.init();

        $(list)
            .on('etalia.list.load', function(e, data) {
                if (data.hasOwnProperty('controlsStates')) {
                    controls.setStates(data.controlsStates);
                }
            });

        $body
            .on(
                'etalia.control.search.change ' +
                'etalia.control.timespan.change ' +
                'etalia.control.cluster.change ' +
                'etalia.control.pinned.change ' +
                'etalia.control.filters.change',
            function() {
                list.load(controls.getStates());
            })
            .on('etalia.publication.pin', function(e, result) {
                $('.thumb[data-paper-id=' + result.getId() + ']')
                    .find('.thumb-pin')
                    .toggleClass('active', result.isPinned());
            })
            .on('etalia.publication.ban', function(e, result) {
                if (result.isBanned()) {
                    $('.thumb[data-paper-id=' + result.getId() + ']').remove();
                }
            })
            .on('etalia.publication.trash', function(e, result) {
                if (result.isTrashed()) {
                    $('.thumb[data-paper-id=' + result.getId() + ']').remove();
                }
            })
            .on('etalia.detail.loading', function() {
                layout.setBusy();
            })
            .on('etalia.detail.loaded', function() {
                layout.setAvailable();
            });

        $('#list').on('click', '.thumb .title a', function(e) {
            detail.load($(e.target).closest('.thumb'));

            e.preventDefault();
            return false;
        });

        // Endless scroll
        $.endlessPaginate({
            paginateOnScroll: true,
            paginateOnScrollMargin: 90,
            onClick: function (context) {
                context.extraData = JSON.stringify(controls.getStates());
                context.urlOrigin = window.location.href;
            },
            onCompleted: function () {
                // Build altmetric badges
                _altmetric_embed_init();
            }
        });
    });
});
