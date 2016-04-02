define([
    'jquery',
    'app/ui/layout',
    'app/ui/controls',
    'bootstrap'
], function($, layout, controls) {

    var List = function(options) {

        this.config = $.extend({
            debug: false,
            element: '#list',
            container: '.thumb-list'
        }, options);

        this.$element = $(this.config.element);

        this.loadXhr = null;
    };

});
