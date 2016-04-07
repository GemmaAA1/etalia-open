define([
    'app',
    'text!app/templates/thread/detail.html',
    'app/view/user/thumb',
    'app/view/user/list',
    'app/view/thread/post/list',
    'app/view/thread/form-edit'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.Detail = App.Backbone.View.extend({

        tagName: 'div',
        className: 'inner',

        template: App.Handlebars.compile(template),

        events: {
            "click .thread-edit": "onEditClick"
        },

        initialize: function () {
            //this.listenTo(this.model, "change:state", this.onThreadStateChange);
            this.listenTo(this.model, "sync", this.render);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var $contents = this.$('.thread-content, .thread-info').hide();

            var form = App.View.Thread.EditForm.create({
                model: this.model
            },{
                $target: this.$('[data-thread-edit-form]')
            });

            this.listenToOnce(form, 'validation_success', function() {
                form.model.save();
            });
            this.listenToOnce(form, 'cancel', function() {
                form.$el.replaceWith('<div data-thread-edit-form></div>');
                form.remove();
                $contents.show();
            });
        },

        render: function() {
            App.log('ThreadDetailView::render');

            this.$el.html(this.template(this.model.attributes));

            // User thumb
            App.View.User.Thumb.create({
                model: this.model.get('user')
            }, {
                $target: this.$('[data-user-placeholder]')
            });

            // Members list
            App.View.User.List.create({
                model: this.model.get('members')
            }, {
                $target: this.$('[data-members-placeholder]')
            });

            // Posts list
            App.View.Thread.PostList.create({
                thread: this.model
            }, {
                $target: this.$('[data-posts-placeholder]')
            });

            return this;
        }
    });
});
