define([
    'jquery',
    'app',
    'text!app/templates/detail.html',
    'app/model/detail'
], function ($, App, template) {

    var defaults = {};

    return App.View.Detail = App.Backbone.View.extend({
        tagName: 'div',
        id: 'detail',

        template: App.Handlebars.compile(template),

        events: {
            "click #detail-close": "close",
            "click #detail-prev": "prev",
            "click #detail-next": "next"
        },

        initialize: function (options) {
            App.defaults(options, defaults);

            if (!this.model instanceof App.Model.Detail) {
                throw 'Unexpected model type';
            }
        },

        render: function () {
            var prev = this.model.get('prev'),
                next = this.model.get('next');

            $('detail-placeholder').replaceWith(this.$el.html(this.template({
                prev: prev ? prev.attributes : null,
                next: next ? next.attributes : null
            })));

            this.model
                .get('view').render().$el
                .appendTo(this.$('.document .wrapper'));

            return this.open();
        },

        prev: function() {
            var prev = this.model.get('prev');
            if (prev) {
                this.trigger('detail:prev', prev, this);
            }
        },

        next: function() {
            var next = this.model.get('next');
            if (next) {
                this.trigger('detail:next', next, this);
            }
        },

        open: function() {
            $('body').addClass('detail-opened');

            return this;
        },

        close: function() {
            $('body').removeClass('detail-opened');

            return this;
        }
    });
});
